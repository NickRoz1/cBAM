module reader;

import std.stdio;
import std.exception;
import std.bitmanip;
import std.conv;
import std.traits;
import std.algorithm.iteration;
import std.range.primitives;
import std.file;
import std.algorithm.searching;
import std.datetime.stopwatch;

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
            parseColumnChunk(recordBuf, columnType, plainColumn);
        }

        return recordBuf;
    }

    void parseColumnChunk(RawReadBlob[] recordBuf, ColumnTypes columnType, ubyte[] plainColumn){
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

    enum bamFields{refID, pos, bin, mapq, flag, nextRefID,
        nextPos, tlen, read_name, cigar, raw_qual, raw_sequence}

    T getColumn(T)(bamFields bamField){
        uint num_rows = reduce!((a,b) => a + b)(map!(row_group => row_group.num_rows)
                            (fileMeta.rowGroups));
        
        static if(is(T==int[])){
            switch(bamField){
                case bamFields.refID:{
                    byte[] buf = getRawColumn(ColumnTypes._refID);
                    T column;
                    column.reserve(num_rows);
                    while(buf.length != 0){
                        column ~= read!(ElementType!T, Endian.bigEndian, byte[])(buf);
                    }
                    return column;
                }
                case bamFields.pos:{
                    byte[] buf = getRawColumn(ColumnTypes._pos);
                    T column;
                    column.reserve(num_rows);
                    while(buf.length != 0){
                        column ~= read!(ElementType!T, Endian.bigEndian, byte[])(buf);
                    }
                    writeln(column[0]);
                    return column;
                }
                case bamFields.nextRefID:{
                    byte[] buf = getRawColumn(ColumnTypes._next_refID);
                    T column;
                    column.reserve(num_rows);
                    while(buf.length != 0){
                        column ~= read!(ElementType!T, Endian.littleEndian, byte[])(buf);
                    }
                    return column;
                }
                case bamFields.nextPos:{
                    byte[] buf = getRawColumn(ColumnTypes._next_pos);
                    T column;
                    column.reserve(num_rows);
                    while(buf.length != 0){
                        column ~= read!(ElementType!T, Endian.littleEndian, byte[])(buf);
                    }
                    return column;
                }
                case bamFields.tlen:{
                    byte[] buf = getRawColumn(ColumnTypes._tlen);
                    T column;
                    column.reserve(num_rows);
                    while(buf.length != 0){
                        column ~= read!(ElementType!T, Endian.littleEndian, byte[])(buf);
                    }
                    return column;
                }
                default:{
                    assert(false, "This field is hold in other type");
                }
            }
        }
        else static if(is(T==short[])){
            switch(bamField){
                case bamFields.bin:{
                    byte[] buf = getRawColumn(ColumnTypes._bin_mq_nl);
                    T column;
                    column.reserve(num_rows);
                    while(buf.length != 0){
                        buf = buf[ushort.sizeof..$];
                        column ~= read!(ElementType!T, Endian.littleEndian, byte[])(buf);
                    }
                    return column;
                }
                case bamFields.flag:{
                    byte[] buf = getRawColumn(ColumnTypes._flag_nc);
                    T column;
                    column.reserve(num_rows);
                    while(buf.length != 0){
                        buf = buf[ushort.sizeof..$];
                        column ~= read!(ElementType!T, Endian.littleEndian, byte[])(buf);
                    }
                    return column;
                }
                default:{
                    assert(false, "This field is hold in other type");
                }
            }
        }
        else static if(is(T==ubyte[])){
            assert(bamField == bamFields.mapq, "This field is hold in other type");
            byte[] buf = getRawColumn(ColumnTypes._bin_mq_nl);
            T column;
            column.reserve(num_rows);
            while(buf.length != 0){
                buf = buf[ubyte.sizeof..$];
                column ~= read!(ElementType!T, Endian.littleEndian, byte[])(buf);
                buf = buf[ushort.sizeof..$];
            }
            return column;
        }
        else static if(is(T==ubyte[][])){
            switch(bamField){
                case bamFields.read_name:{
                    byte[] buf = getRawColumn(ColumnTypes.read_name);
                    T column;
                    column.reserve(num_rows);
                    int offset = 0;
                    while(buf.length != 0){
                        int length = 
                            peek!(int, Endian.littleEndian, byte[])(buf, &offset);
                        column ~= buf[offset..offset+length];
                        offset += length;
                    }
                    return column;
                }
                case bamFields.cigar:{
                    byte[] buf = getRawColumn(ColumnTypes.raw_cigar);
                    T column;
                    column.reserve(num_rows);
                    int offset = 0;
                    while(buf.length != 0){
                        int length = 
                            peek!(int, Endian.littleEndian, byte[])(buf, &offset);
                        column ~= buf[offset..offset+length];
                        offset += length;
                    }
                    return column;
                }
                case bamFields.raw_qual:{
                    byte[] buf = getRawColumn(ColumnTypes.raw_qual);
                    T column;
                    column.reserve(num_rows);
                    int offset = 0;
                    while(buf.length != 0){
                        int length = 
                            peek!(int, Endian.littleEndian, byte[])(buf, &offset);
                        column ~= buf[offset..offset+length];
                        offset += length;
                    }
                    return column;
                }
                case bamFields.raw_sequence:{
                    byte[] buf = getRawColumn(ColumnTypes.raw_sequence);
                    T column;
                    column.reserve(num_rows);
                    int offset = 0;
                    while(buf.length != 0){
                        int length = 
                            peek!(int, Endian.littleEndian, byte[])(buf, &offset);
                        column ~= buf[offset..offset+length];
                        offset += length;
                    }
                    return column;
                }
                default:{
                    assert(false, "This field is hold in other type");
                }
            }
        }

    }

    
    
    byte[] getRawColumn(ColumnTypes columnType){
        ulong raw_column_size = reduce!((a,b) => a + b)
                (map!(row_group => row_group.columnsSizes[columnType])
                    (fileMeta.rowGroups));

        ulong max_column_chunk = maxElement(map!(row_group => row_group.columnsSizes[columnType])
                    (fileMeta.rowGroups));
        byte[] buffer;
        buffer.length = raw_column_size;
        //buffer.reserve(max_column_chunk);

        
        // Dangerous, may produce wrong results. Investigate
        // auto ok = snappy_uncompressed_length(buffer.ptr, raw_column_size, &prediction);
        // enforce(ok == snappy_status.SNAPPY_OK);
        // byte[] outBuf;
        // writeln(prediction);

        ulong offset = 0;
        foreach(rowGroup; fileMeta.rowGroups){
            file.seek(rowGroup.columnsOffsets[columnType]);
            ulong size = rowGroup.columnsSizes[columnType];
            file.rawRead(buffer[offset..offset+size]);
            offset += size;
            //outBuf ~= Snappy.uncompress(buffer);
        }
        size_t prediction;
        ulong offset1 = 0;
        foreach(rowGroup; fileMeta.rowGroups){
            size_t temp;
            auto tempSlice = buffer[offset1..offset1+rowGroup.columnsSizes[columnType]];
            snappy_uncompressed_length(tempSlice.ptr, tempSlice.length, &temp);
            prediction += temp;
        }
        //writeln(prediction);
        //auto ok = snappy_uncompressed_length(buffer.ptr, raw_column_size, &prediction);
        //enforce(ok == snappy_status.SNAPPY_OK);
        byte[] outBuf;
        //writeln(prediction);
        outBuf.reserve(prediction);

        foreach(rowGroup; fileMeta.rowGroups){
            outBuf ~= Snappy.uncompress(buffer[0..rowGroup.columnsSizes[columnType]]);
            buffer = buffer[rowGroup.columnsSizes[columnType]..$];
        }

        //writeln("TEST");
        //writeln(outBuf.length);

        return outBuf;
    }

    unittest{
        auto fn = getcwd() ~ "/source/tests/test3.cbam";
        File file = File(fn, "r");
        FileReader fileR = new FileReader(file);
        auto sw = StopWatch(AutoStart.no);
        sw.start();
        int[] test = fileR.getColumn!(int[])(bamFields.pos);
        sw.stop();
        writeln(test.length);

        //test.length = 0;
        writeln("CBAM DURATION:");
        writeln(sw.peek.total!"usecs");
        fn = getcwd() ~ "/source/tests/test2.bam";
        auto reader = BamBlobReader(fn);
        int i = 0;
        sw.start();
        while(!reader.empty){
            test[i] = reader.fetch().pos;
            ++i;
        }
        sw.stop();
        writeln("BAM DURATION:");
        writeln(sw.peek.total!"usecs");
        writeln(test.length);
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