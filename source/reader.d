module reader;

import std.stdio;
import std.exception;
import std.bitmanip;
import std.conv;
import std.traits;

import source.writer;
import bio.std.experimental.hts.bam.header;
import utils.bam.reader;
import snappy.snappy;



class FileReader {
    File file;
    FileMeta fileMeta;
    BamHeader bamHeader;

    this(File f) {
        file = f;
        parse_meta();
    }

    ~this(){
        if(file.isOpen) return;
        file.close();
    }

    RawReadBlob[] readRowGroup(int i){
        enforce(i >= 0 && i < fileMeta.rowGroups.length, "No such rowgroup " ~ to!string(i));

        RowGroupMeta rowGroup = fileMeta.rowGroups[i];
        RawReadBlob[] recordBuf;
        recordBuf.length = rowGroup.num_rows;
        
        file.seek(rowGroup.columnsOffsets[0]);
        ubyte[] buf;
        buf.length = rowGroup.total_byte_size;
        file.rawRead(buf);

        
        ulong offset = 0;
        foreach(columnType; EnumMembers!ColumnTypes){
            ulong compressedSize = rowGroup.columnsSizes[columnType];
            ubyte[] plainColumn = 
                cast(ubyte[])Snappy.uncompress(cast(byte[])buf[offset..offset + compressedSize]);
            offset += compressedSize;
            
            // _blob_size should come before raw fields
            parseColumn(recordBuf, columnType, plainColumn);
        }

        return recordBuf;
    }

    void parseColumn(RawReadBlob[] recordBuf, ColumnTypes columnType, ubyte[] plainColumn){
        pragma(inline, true);

        void injectField(RawReadBlob readBlob, int offset, int size, ubyte[] data){
            assert(readBlob._data.length >= offset + size, "_blob is smaller than raw data");
            readBlob._data[offset..offset+size] = data.dup; // Avoid dup? Ensure buf is immutable then
        }

        switch(columnType){
                case ColumnTypes._refID: {
                    foreach(ref record; recordBuf){
                        record.refid = read!(int, Endian.bigEndian, ubyte[])
                            (plainColumn);
                    }
                    break;
                }
                case ColumnTypes._pos: {
                    foreach(ref record; recordBuf){
                        record.pos = read!(uint, Endian.bigEndian, ubyte[])
                            (plainColumn);
                    }
                    break;
                }
                case ColumnTypes._blob_size: {
                    foreach(ref record; recordBuf){
                        record._data.length = read!(uint, Endian.bigEndian, ubyte[])
                            (plainColumn);
                    }
                    break;
                }
                case ColumnTypes._bin_mq_nl: {
                    ulong columnOffset = 0;
                    ubyte[] data;
                    foreach(ref record; recordBuf){
                        data = cast(ubyte[])plainColumn[columnOffset..columnOffset+simple_field_size];
                        injectField(record, Offset.bin_mq_nl, simple_field_size, data);
                        columnOffset += simple_field_size;
                    }
                    break;
                }
                case ColumnTypes._flag_nc: {
                    ulong columnOffset = 0;
                    ubyte[] data;
                    foreach(ref record; recordBuf){
                        data = cast(ubyte[])plainColumn[columnOffset..columnOffset+simple_field_size];
                        injectField(record, Offset.flag_nc, simple_field_size, data);
                        columnOffset += simple_field_size;
                    }
                    break;
                }
                case ColumnTypes.sequence_length: {
                    ulong columnOffset = 0;
                    ubyte[] data;
                    foreach(ref record; recordBuf){
                        data = cast(ubyte[])plainColumn[columnOffset..columnOffset+simple_field_size];
                        injectField(record, Offset.l_seq, simple_field_size, data);
                        columnOffset += simple_field_size;
                    }
                    break;
                }
                case ColumnTypes._next_refID: {
                    ulong columnOffset = 0;
                    ubyte[] data;
                    foreach(ref record; recordBuf){
                        data = cast(ubyte[])plainColumn[columnOffset..columnOffset+simple_field_size];
                        injectField(record, Offset.next_refID, simple_field_size, data);
                        columnOffset += simple_field_size;
                    }
                    break;
                }
                case ColumnTypes._next_pos: {
                    ulong columnOffset = 0;
                    ubyte[] data;
                    foreach(ref record; recordBuf){
                        data = cast(ubyte[])plainColumn[columnOffset..columnOffset+simple_field_size];
                        injectField(record, Offset.next_pos, simple_field_size, data);
                        columnOffset += simple_field_size;
                    }
                    break;
                }
                case ColumnTypes._tlen: {
                    ulong columnOffset = 0;
                    ubyte[] data;
                    foreach(ref record; recordBuf){
                        data = cast(ubyte[])plainColumn[columnOffset..columnOffset+simple_field_size];
                        injectField(record, Offset.tlen, simple_field_size, data);
                        columnOffset += simple_field_size;
                    }
                    break;
                }
                case ColumnTypes.read_name: {
                    size_t columnOffset = 0;
                    ubyte[] data;
                    foreach(ref record; recordBuf){
                        int len = peek!(int, Endian.littleEndian, ubyte[])
                            (plainColumn, &columnOffset);
                        data = cast(ubyte[])plainColumn[columnOffset..columnOffset+len];
                        injectField(record, Offset.read_name, len, data);
                        columnOffset += len;
                    }
                    break;
                }
                case ColumnTypes.raw_cigar: {
                    ulong columnOffset = 0;
                    ubyte[] data;
                    foreach(ref record; recordBuf){
                        int len = peek!(int, Endian.littleEndian, ubyte[])
                            (plainColumn, &columnOffset);
                        data = cast(ubyte[])plainColumn[columnOffset..columnOffset+len];
                        injectField(record, record._cigar_offset, len, data);
                        columnOffset += len;
                    }
                    break;
                }
                case ColumnTypes.raw_qual: {
                    ulong columnOffset = 0;
                    ubyte[] data;
                    foreach(ref record; recordBuf){
                        int len = peek!(int, Endian.littleEndian, ubyte[])
                            (plainColumn, &columnOffset);
                        data = cast(ubyte[])plainColumn[columnOffset..columnOffset+len];
                        injectField(record, record._qual_offset, len, data);
                        columnOffset += len;
                    }
                    break;
                }
                case ColumnTypes.raw_sequence: {
                    ulong columnOffset = 0;
                    ubyte[] data;
                    foreach(ref record; recordBuf){
                        int len = peek!(int, Endian.littleEndian, ubyte[])
                            (plainColumn, &columnOffset);
                        data = cast(ubyte[])plainColumn[columnOffset..columnOffset+len];
                        injectField(record, record._seq_offset, len, data);
                        columnOffset += len;
                    }
                    break;
                }
                default: assert(false, "No such type exists");
            }
    }





    /// Parse FileMeta
    void parse_meta() {
        ulong file_size = file.size();
        enforce(file_size >= CBAM_MAGIC.sizeof, "File size is smaller than MAGIC string");

        file.seek(-1*(2*uint.sizeof+ulong.sizeof), SEEK_END);
        ubyte[] buffer;
        buffer.length = 16;
        file.rawRead(buffer);
        
        
        enforce(buffer[$-uint.sizeof..$] == CBAM_MAGIC, "Invalid file");

        ulong meta_offset = read!(ulong, Endian.littleEndian, ubyte[])
                                (buffer);
        uint meta_size = read!(uint, Endian.littleEndian, ubyte[])
                                (buffer);
        
        file.seek(meta_offset);
        if(meta_size > 0){
            buffer.length = meta_size;
            file.rawRead(buffer);

            uint rowGroupsAmount = read!(int, Endian.littleEndian, ubyte[])(buffer);

            fileMeta.rowGroups.length = rowGroupsAmount;

            foreach(ref rowGroup; fileMeta.rowGroups){
                foreach(ref columnOffset; rowGroup.columnsOffsets){
                    columnOffset = read!(ulong, Endian.bigEndian, ubyte[])(buffer);
                }

                foreach(ref columnSize; rowGroup.columnsSizes){
                    columnSize = read!(ulong, Endian.bigEndian, ubyte[])(buffer);
                }

                rowGroup.total_byte_size = read!(ulong, Endian.bigEndian, ubyte[])(buffer);

                rowGroup.num_rows = read!(uint, Endian.bigEndian, ubyte[])(buffer);

            }
        }
        parseBamHeader();
    }

    void parseBamHeader() {
        ubyte[] buffer;
        buffer.length = uint.sizeof;
        file.rawRead(buffer);
        uint header_text_size = read!(uint, Endian.littleEndian, ubyte[])(buffer);

        if(header_text_size > 0){
            buffer.length = header_text_size;
            auto test = file.rawRead(buffer);
            bamHeader.text = cast(string)buffer.dup;
        }
    
        buffer.length = uint.sizeof;
        file.rawRead(buffer);
        uint num_ref_seq = read!(uint, Endian.littleEndian, ubyte[])(buffer);
        bamHeader.refs.length = num_ref_seq;

        foreach(ref ref_seq; bamHeader.refs){
            buffer.length = uint.sizeof;

            file.rawRead(buffer);
            uint name_length = read!(uint, Endian.littleEndian, ubyte[])(buffer);
            if(name_length > 0){
                buffer.length = name_length;
                file.rawRead(buffer);
                ref_seq.name = cast(string)buffer.dup;
            }

            buffer.length = ulong.sizeof;
            file.rawRead(buffer);
            ref_seq.length = read!(ulong, Endian.littleEndian, ubyte[])(buffer);
        }
    }
    
}