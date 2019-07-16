module source.writer;
import std.stdio;
import std.string;
import std.traits;
import std.conv;
import std.file;
import std.bitmanip;
import std.algorithm;
import std.algorithm.iteration;
import std.exception;

//import sambamba.utils.lz4;


import utils.bam.reader;
import bio.std.experimental.hts.bam.header;
import snappy.snappy;
import reader;

const ubyte[4] CBAM_MAGIC = ['C','B','A','M'];
const int BATCH_SIZE = 1_000_000;
const uint simple_field_size = uint.sizeof;

int bamToCbam(const string fn, const string ofn){
    //File outfile = File(ofn);
    BamReadBlobStream reader = BamReadBlobStream(fn); 
    File file = File(ofn, "w");

    FileWriter fileWriter = new FileWriter(file, reader.header);
    int recordCount = 0;
    RawReadBlob[] recordBuf;
    recordBuf.length = BATCH_SIZE;
    while(!reader.empty()){
        uint num_rows = 0;
        
        while(num_rows < BATCH_SIZE && !reader.empty()){
            recordBuf[num_rows] = reader.front();
            reader.popFront();
            ++num_rows;
        }

        if(reader.empty()){ // 
            recordBuf[num_rows] = reader.front();
            ++num_rows;
        }

        recordCount += num_rows;
        
        fileWriter.writeRowGroup(recordBuf, num_rows);
    }
    fileWriter.close();

    return recordCount;
}

struct RowGroupMeta{
    ulong[EnumMembers!ColumnTypes.length] columnsOffsets;
    ulong[EnumMembers!ColumnTypes.length] columnsSizes;
    ulong total_byte_size;
    uint num_rows;
}

struct FileMeta{
    RowGroupMeta[] rowGroups;
    //uint compressionLevel;
}

// Order matters
enum ColumnTypes {_refID, _pos, _blob_size, _bin_mq_nl, _flag_nc, sequence_length, _next_refID, _next_pos,
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
                    writeFieldToBuf(buf, columnType, recordBuf[i], i*simple_field_size);
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

    unittest {
        auto fn = getcwd() ~ "/source/tests/test1.bam";

        File testFile = File.tmpfile();
        scope(exit) testFile.close();
        
        BamBlobReader reader = BamBlobReader(fn); 

        FileWriter fileWriter = new FileWriter(testFile, reader.header);

        const uint test_batch_size = 10;
        RawReadBlob[] recordBuf;
        recordBuf.length = BATCH_SIZE;
        uint num_rows = 0;
        while(!reader.empty()){
            num_rows = 0;
            
            while(num_rows < BATCH_SIZE && !reader.empty()){
                recordBuf[num_rows] = reader.fetch();
                ++num_rows;
            }
            
            fileWriter.writeRowGroup(recordBuf, num_rows);
        }
        fileWriter.writeMeta();
        testFile.flush();
        testFile.rewind();
        
        FileReader fileReader = new FileReader(testFile);
        BamBlobReader reader2 = BamBlobReader(fn);
        int rowGroupNum = 0;
        while(!reader.empty()){
            while(num_rows < BATCH_SIZE && !reader.empty()){
                recordBuf[num_rows] = reader.fetch();
                ++num_rows;
            }
            auto testBuf = fileReader.readRowGroup(rowGroupNum);
            for(int i = 0; i < num_rows; i++){
                assert(recordBuf[i]._bin_mq_nl == testBuf[i]._bin_mq_nl);
                assert(recordBuf[i]._flag_nc == testBuf[i]._flag_nc);
                assert(recordBuf[i].sequence_length == testBuf[i].sequence_length);
                assert(recordBuf[i]._next_pos == testBuf[i]._next_pos);
                assert(recordBuf[i]._next_refID == testBuf[i]._next_refID);
                assert(recordBuf[i]._mapq == testBuf[i]._mapq);
                assert(recordBuf[i]._flag == testBuf[i]._flag);
                assert(equal(recordBuf[i].read_name, testBuf[i].read_name));
                assert(equal(recordBuf[i].raw_cigar, testBuf[i].raw_cigar));
                assert(equal(recordBuf[i].raw_qual, testBuf[i].raw_qual));
                assert(equal(recordBuf[i].raw_sequence, testBuf[i].raw_sequence));
            }
            rowGroupNum++;
        }
    }

    void writeFieldToBuf(ubyte[] buf, ColumnTypes columnType, RawReadBlob readBlob, int offset){
        //pragma(inline, true);


        switch(columnType) { 
            case ColumnTypes._refID: {
                std.bitmanip.write!(int, Endian.littleEndian, ubyte[])
                    (buf, readBlob.refid, offset);
                break;
            }
            case ColumnTypes._pos: {
                std.bitmanip.write!(int, Endian.littleEndian, ubyte[])
                    (buf, readBlob.pos, offset);
                break;
            }
            case ColumnTypes._blob_size: {
                uint blob_size = cast(int)readBlob._data.length;
                std.bitmanip.write(buf, blob_size, offset);
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
        //pragma(inline, true);
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
        //pragma(inline, true);    

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
        //pragma(inline, true);

        byte[] buf = Snappy.compress(cast(byte[])column);
        //auto buf = column;
        file.rawWrite(buf);
        return buf.length;
    }

    //       ...|fileMeta|bamHeader| meta_offset |meta_size| MAGIC |
    //       ...|        |         |   8 bytes   | 4 bytes |4 bytes|
    //       ...|        |         |             |         |       |
    
    /// Write meta to file
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

        writeBamHeader(bamHeader, file);

        buf.length = ulong.sizeof;
        writeToFile(meta_offset, file, buf);
        writeToFile(meta_size, file, buf);
        file.rawWrite(CBAM_MAGIC);
    }

    unittest{
        auto fn = getcwd() ~ "/source/tests/test1.bam";
        
        File testFile = File.tmpfile();
        scope(exit) testFile.close();
        
        BamBlobReader reader = BamBlobReader(fn); 

        FileWriter fileWriter = new FileWriter(testFile, reader.header);
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

        fileWriter.writeMeta();
        testFile.flush();
        
        testFile.rewind();
        FileReader fileReader = new FileReader(testFile);

        assert(fileReader.fileMeta == fileWriter.fileMeta);
        assert(equal!"a == b"(fileReader.bamHeader.refs, fileWriter.bamHeader.refs)); 
        assert(equal(fileReader.bamHeader.text,fileWriter.bamHeader.text)); 
    }

    /// Writes BamHeader to the file
    static void writeBamHeader(BamHeader bamHeader, File file){
        //pragma(inline, true);

        ubyte[] buf;
        buf.length = ulong.sizeof;

        auto text_size = cast(uint)bamHeader.text.length;

        writeToFile(text_size, file, buf);
        file.rawWrite(cast(ubyte[])bamHeader.text);
        writeToFile(cast(uint)bamHeader.refs.length, file, buf);

        foreach(ref_seq; bamHeader.refs) {
            ubyte[] ref_name = cast(ubyte[])ref_seq.name.dup;
            auto ref_name_length = cast(uint)ref_name.length;

            writeToFile(ref_name_length, file, buf);
            file.rawWrite(ref_name);
            writeToFile(ref_seq.length, file, buf);
        }
    }

    unittest{
        File testFile = File.tmpfile();
        scope(exit) testFile.close();

        auto fn = getcwd() ~ "/source/tests/test1.bam";
        auto bamReader = BamBlobReader(fn);
        BamHeader bamHeader = bamReader.header;
        
        writeBamHeader(bamHeader, testFile);

        ubyte[] buf;
        buf.length = ulong.sizeof;
        writeToFile(cast(ulong)0, testFile, buf);
        writeToFile(cast(uint)0, testFile, buf);
        testFile.rawWrite(CBAM_MAGIC);
        testFile.flush();
        
        testFile.rewind();
        auto fileReader = new FileReader(testFile);

        // writeln(bamHeader.refs.ptr);
        // writeln(fileReader.bamHeader.refs.ptr);
        assert(equal!"a == b"(fileReader.bamHeader.refs, bamHeader.refs), "Refs are corrupted");
        assert(equal(fileReader.bamHeader.text, bamHeader.text), "Text is corrupted");
    }
    
    /// Write T obj to file. buf.length must be <= ulong.sizeof
    static void writeToFile(T)(T obj, File file, ubyte[] buf){
        enforce(T.sizeof <= ulong.sizeof);
        enforce(buf.length <= ulong.sizeof);
        enforce(T.sizeof <= buf.length, "No types longer than ulong"); 
        
        write!(T, Endian.littleEndian, ubyte[])
                (buf, obj, 0);
        file.rawWrite(buf[0..T.sizeof]);
    }
    unittest{
        File testFile = File.tmpfile();
        scope(exit) testFile.close();


        ubyte[] buf;
        buf.length = ulong.sizeof;

        const ubyte ub = 10;
        const ushort us = 12;
        const uint ui = 14;
        const ulong ul = 16;

        writeToFile(ub, testFile, buf);
        writeToFile(us, testFile, buf);
        writeToFile(ui, testFile, buf);
        writeToFile(ul, testFile, buf);
        testFile.flush();

        testFile.rewind();
        buf.length = ubyte.sizeof + ushort.sizeof + uint.sizeof + ulong.sizeof;
        testFile.rawRead(buf);

        assert(read!(ubyte, Endian.littleEndian, ubyte[])(buf) == ub);
        assert(read!(ushort, Endian.littleEndian, ubyte[])(buf) == us);
        assert(read!(uint, Endian.littleEndian, ubyte[])(buf) == ui);
        assert(read!(ulong, Endian.littleEndian, ubyte[])(buf) == ul);
    }
}


