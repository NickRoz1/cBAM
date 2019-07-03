module source.writer;
import std.stdio;
import std.string;
import std.traits;
import std.bitmanip;
import std.algorithm.iteration;
import std.exception;

//import sambamba.utils.lz4;


import utils.bam.reader;
import bio.std.experimental.hts.bam.header;
import snappy.snappy;

const ubyte[4] CBAM_MAGIC = ['C','B','A','M'];
const int BATCH_SIZE = 100_000;
const uint simple_field_size = uint.sizeof;

int bamToCbam(const string fn, const string ofn){
    //File outfile = File(ofn);
    BamBlobReader reader = BamBlobReader(fn); 
    File file = File(ofn, "w");

    FileWriter fileWriter = new FileWriter(file, reader.header);
    int recordCount = 0;
    RawReadBlob[] recordBuf;
    recordBuf.length = BATCH_SIZE;
    while(!reader.empty()){
        uint num_rows = 0;
        
        while(num_rows < BATCH_SIZE && !reader.empty()){
            recordBuf[num_rows] = reader.fetch();
            ++num_rows;
        }

        recordCount += num_rows;
        
        fileWriter.writeRowGroup(recordBuf, num_rows);
    }
    fileWriter.close();

    return recordCount;
}

struct RowGroupMeta{
    ulong[12] columnsOffsets;
    ulong[12] columnsSizes;
    ulong total_byte_size;
    uint num_rows;
}

struct FileMeta{
    RowGroupMeta[] rowGroups;
    //uint compressionLevel;
}

// Order matters - raw qual buffer reused for raw sequence
enum ColumnTypes {_refID, _pos, _bin_mq_nl, _flag_nc, sequence_length, _next_refID, _next_pos,
_tlen, read_name, raw_cigar, raw_qual, raw_sequence}

class FileWriter{
    File file;
    FileMeta fileMeta;
    BamHeader bamHeader;
    uint curOffset;
    
    this(File fn, BamHeader bamHead){
        file = fn;
        file.rawWrite(CBAM_MAGIC);
        curOffset = 0;
        bamHeader = bamHead;         
    }
    
    ~this() {
        close();
    }

    void close(){
        if(!file.isOpen) return;
        writeMeta();
        file.close(); 
    }
    void writeRowGroup(RawReadBlob[] recordBuf, uint num_rows){
        assert(file.isOpen());
        
        RowGroupMeta rowGroupMeta;
        rowGroupMeta.num_rows = num_rows;

        uint total_size = 0;
        ubyte[] buf;
        buf.length = num_rows * int.sizeof; // biggest field is int32

        foreach(columnType; EnumMembers!ColumnTypes) {
            if(columnType < ColumnTypes.read_name) { // Types which memory size is known
                rowGroupMeta.columnsOffsets[columnType] = file.tell; // may be wrong, check
                for(int i = 0; i < num_rows; ++i) {
                    writeFieldToBuf(buf, columnType, recordBuf[i], i);
                }
                rowGroupMeta.columnsSizes[columnType] = writeColumn(buf[0..num_rows*simple_field_size]);
            }
            else {
                if(columnType != ColumnTypes.raw_sequence) { // raw_qual >= raw_sequence
                    buf.length = calcBufSize(columnType, recordBuf)+int.sizeof*num_rows;
                }
                rowGroupMeta.columnsOffsets[columnType] = file.tell();

                int currentPos = 0;
                for(int i = 0; i < num_rows; ++i) {
                    writeVarFieldToBuf(buf, columnType, recordBuf[i], currentPos);
                }
                
                rowGroupMeta.columnsSizes[columnType] = writeColumn(buf[0..currentPos]);
            }
            //writeln(columnType);
        }

        rowGroupMeta.total_byte_size = reduce!((a,b) => a + b)(rowGroupMeta.columnsSizes);
        fileMeta.rowGroups ~= rowGroupMeta;
    }

    void writeFieldToBuf(ubyte[] buf, ColumnTypes columnType, RawReadBlob readBlob, int offset){
        pragma(inline, true);


        switch(columnType) { 
            case ColumnTypes._refID: {
                std.bitmanip.write(buf, readBlob.refid, offset);
                break;
            }
            case ColumnTypes._pos: {
                std.bitmanip.write(buf, readBlob.pos, offset);
                break;
            }
            case ColumnTypes._bin_mq_nl: {
                buf[offset..offset+simple_field_size] = readBlob.raw_bin_mq_nl;
                break;
            }
            case ColumnTypes.sequence_length: {
                buf[offset..offset+simple_field_size] = readBlob.raw_sequence_length;
                break;
            }
            case ColumnTypes._flag_nc: {
                buf[offset..offset+simple_field_size] = readBlob.raw_flag_nc;
                break;
            }
            case ColumnTypes._next_pos: {
                buf[offset..offset+simple_field_size] = readBlob.raw_next_pos;
                break;
            }
            case ColumnTypes._next_refID: {
                buf[offset..offset+simple_field_size] = readBlob.raw_next_refID;
                break;
            }
            case ColumnTypes._tlen: {
                buf[offset..offset+simple_field_size] = readBlob.raw_tlen;
                break;
            }
            default:{ assert(false, "No such type exists"); } //should throw
        }
    }

    void writeVarFieldToBuf(ubyte[] buf, ColumnTypes columnType, RawReadBlob readBlob,
                            ref int offset){
        pragma(inline, true);
        //if(!(offset%10000)) writeln(offset);

        switch(columnType) { 
            case ColumnTypes.read_name: {
                int fieldSize = cast(int)readBlob.read_name.length;
                write!(int, Endian.littleEndian, ubyte[])(buf, fieldSize, offset);
                offset += int.sizeof;
                buf[offset..offset + readBlob.read_name.length] = readBlob.read_name;
                offset += readBlob.read_name.length;
                break;
            }
            case ColumnTypes.raw_cigar: {
                int fieldSize = cast(int)readBlob.raw_cigar.length;
                write!(int, Endian.littleEndian, ubyte[])(buf, fieldSize, offset);
                offset += int.sizeof;
                buf[offset..offset + readBlob.raw_cigar.length] = readBlob.raw_cigar;
                offset += readBlob.raw_cigar.length;
                break;
            }
            case ColumnTypes.raw_qual: {
                int fieldSize = cast(int)readBlob.raw_qual.length;
                write!(int, Endian.littleEndian, ubyte[])(buf, fieldSize, offset);
                offset += int.sizeof;
                buf[offset..offset + readBlob.raw_qual.length] = readBlob.raw_qual;
                offset += readBlob.raw_qual.length;
                break;
            }
            case ColumnTypes.raw_sequence: {
                int fieldSize = cast(int)readBlob.raw_sequence.length;
                write!(int, Endian.littleEndian, ubyte[])(buf, fieldSize, offset);
                offset += int.sizeof;
                buf[offset..offset + readBlob.raw_sequence.length] = readBlob.raw_sequence;
                offset += readBlob.raw_sequence.length;
                break;
            }

            default:{ assert(false, "No such type exists"); }//should throw
        }
    }

    /// Returns size of buffer to allocate to hold columnType field of all records
    uint calcBufSize(ColumnTypes columnType, RawReadBlob[] recordBuf){
        pragma(inline, true);    

        switch(columnType){
            case ColumnTypes.read_name:
                return reduce!((a,b) => a + b)(map!(a => a._l_read_name)(recordBuf));
            case ColumnTypes.raw_cigar:
                return cast(uint)(
                    int.sizeof*reduce!((a,b) => a + b)(map!(a => a._n_cigar_op)(recordBuf)));
            case ColumnTypes.raw_sequence:
                return (reduce!((a,b) => a + b)(map!(a => a.sequence_length)(recordBuf))+1)/2;
            case ColumnTypes.raw_qual:
                return reduce!((a,b) => a + b)(map!(a => a.sequence_length)(recordBuf));
            default: { assert(false); }
        }
    }

    ulong writeColumn(ubyte[] column){
        pragma(inline, true);

        byte[] buf = Snappy.compress(cast(byte[])column);
        //auto buf = column;
        file.rawWrite(buf);
        return buf.length;
    }

    void writeMeta(){
        ulong meta_offset = file.tell();

        ubyte[] buf;
        ulong offset = 0;
        uint rowGroupsAmount = cast(uint)fileMeta.rowGroups.length; 
        buf.length = int.sizeof + rowGroupsAmount * RowGroupMeta.sizeof;
        write!(int, Endian.littleEndian, ubyte[])(buf, rowGroupsAmount, offset);
        offset += int.sizeof;

        foreach(rowGroup; fileMeta.rowGroups) {
            foreach(columnOffset; rowGroup.columnsOffsets){
                std.bitmanip.write(buf, columnOffset, offset);
                offset += ulong.sizeof;
            }
            foreach(columnSize; rowGroup.columnsSizes){
                std.bitmanip.write(buf, columnSize, offset);
                offset += ulong.sizeof;
            }
            
            std.bitmanip.write(buf, rowGroup.total_byte_size, offset);
            offset += ulong.sizeof;

            std.bitmanip.write(buf, rowGroup.num_rows, offset);
            offset += uint.sizeof;
        }

        file.rawWrite(buf);
        uint meta_size = cast(uint)buf.length;

        uint text_length = cast(uint)bamHeader.text.length;
        ubyte[] bufToWrite;
        bufToWrite.length = uint.sizeof;

        write!(uint, Endian.littleEndian, ubyte[])(bufToWrite, text_length, 0);
        file.rawWrite(bufToWrite);
        file.rawWrite(representation(bamHeader.text));

        write!(uint, Endian.littleEndian, ubyte[])
            (bufToWrite, cast(uint)bamHeader.refs.length, 0);
        file.rawWrite(bufToWrite);
        foreach(ref_seq; bamHeader.refs) {
            bufToWrite.length = uint.sizeof; // introduce two buffers for different types
            uint len = cast(uint)ref_seq.name.length;
            write!(uint, Endian.littleEndian, ubyte[])
                (bufToWrite, len, 0);
            file.rawWrite(bufToWrite);
            offset += uint.sizeof;

            file.rawWrite(representation(ref_seq.name));
            offset += len;

            bufToWrite.length = ulong.sizeof;
            write!(ulong, Endian.littleEndian, ubyte[])
                (bufToWrite, ref_seq.length, 0);
            file.rawWrite(bufToWrite);
            offset += ulong.sizeof;
        }

        /*      ...| meta_offset |meta_size| MAGIC |
        /       ...|   8 bytes   | 4 bytes |4 bytes|
        /       ...|             |         |       |
        */

        bufToWrite.length = ulong.sizeof;
        write!(ulong, Endian.littleEndian, ubyte[])
                (bufToWrite, meta_offset, 0);
        file.rawWrite(bufToWrite);

        bufToWrite.length = uint.sizeof;
        write!(uint, Endian.littleEndian, ubyte[])
                (bufToWrite, meta_size, 0);
        file.rawWrite(bufToWrite);
        file.rawWrite(CBAM_MAGIC);
    }
}


