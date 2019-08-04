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
import std.algorithm.comparison;
import std.parallelism;

import source.writer;
//import utils.bam.header;
import bio.std.experimental.hts.bam.header;

//import bio.std.hts.bam.reader;
import utils.bam.reader;
import snappy.snappy;



class FileReader {
    File file;
    FileMeta fileMeta;
    BamHeader bamHeader;

    this(File f) {
        file = f;
        parse_meta();
        //if(fileMeta.rowGroups.length == 0) writeln("File is empty");
    }

    ~this(){
        if(file.isOpen) return;
        file.close();
    }

    void close(){
        file.close();
    }

    // RawReadBlob[] readRowGroup(const int i){
    //     enforce(i >= 0 && i < fileMeta.rowGroups.length, "No such rowgroup " ~ to!string(i));

    //     void injectField(ref RawReadBlob readBlob, int offset, int size, ref ubyte[] data){
    //         assert(readBlob._data.length >= offset + size, "_blob is smaller than raw data");
    //         readBlob._data[offset..offset+size] = data; // dup? Ensure buf is immutable then
    //     }

    //     static int localOff()

    //     void injectField(T) 
    //         if (is(T==Column))
    //         (RawReadBlob[] readsBuf, T col){
    //         static if(is(T == Column!Args, Args...)){
    //             static if(is(Args[0] == ubyte[])){
    //                 alias bamField = Args[1];
    //                 static if(bamField == bamFields.read_name)
    //                     int localOffset = 
    //                 else static if(bamField == bamFields.raw_cigar)

    //                 else static if(bamField == bamFields.raw_sequence)
    //                 else static if(bamField == bamFields.raw_qual)
    //             }
    //             int localOffset = 
    //         }
    //     }
    //     auto sizesCol = getColumn!(int[], bamFields.blob_size).getChunk(i);
        
    //     RawReadBlob[] readsBuf;
    //     readsBuf.length = sizesCol.length;

    //     foreach(i, ref elem; readsBuf){
    //         elem._data.length = sizesCol[i];
    //     }

    //     static foreach(member; EnumMembers!bamFields){
    //         static if(member == bamFields.read_name ||
    //                   member == bamFields.raw_qual  ||
    //                   member == bamFields.raw_sequence ||
    //                   member == bamFields.cigar){
    //             auto col = getColumn!(ubyte[][], member, true).getChunk(i);
    //             injectField(typeof(col))(readsBuf, col);
    //         }
    //         else static if(member != bamFields.mapq && member != bamFields.blob_size){
    //             auto col = getColumn!(int[], member, true).getChunk(i);
    //         }
    //     }
    // }



    enum bamFields{refID, pos, bin, mapq, flag, nextRefID,
        nextPos, tlen, read_name, cigar, raw_sequence, raw_qual, blob_size}

    enum bamfToCbamf = [
        bamFields.refID         : ColumnTypes._refID,
        bamFields.pos           : ColumnTypes._pos,
        bamFields.bin           : ColumnTypes._bin_mq_nl,
        bamFields.mapq          : ColumnTypes._bin_mq_nl,
        bamFields.flag          : ColumnTypes._flag_nc,
        bamFields.nextRefID     : ColumnTypes._next_refID,
        bamFields.nextPos       : ColumnTypes._next_pos,
        bamFields.tlen          : ColumnTypes._tlen,
        bamFields.read_name     : ColumnTypes.read_name,
        bamFields.cigar         : ColumnTypes.raw_cigar,
        bamFields.raw_qual      : ColumnTypes.raw_qual,
        bamFields.raw_sequence  : ColumnTypes.raw_sequence,
        bamFields.blob_size     : ColumnTypes._blob_size
    ];

    Column!(T, bfield) getColumn(T, bamFields bfield)(){
        return Column!(T, bfield)(fileMeta, file);
    }

    Column!(T, bfield, isRaw) getColumn(T, bamFields bfield, bool isRaw)(){
        return Column!(T, bfield, isRaw)(fileMeta, file);
    }

    struct Column(T, bamFields bfield, bool isRaw = false)
        if(is(T == ubyte[][]) || is(T == int[]) 
        || is(T == short[])   || is(T == ubyte[])) {
            
        this(FileMeta meta, File f){
            fileMeta = meta;
            static if(isRaw && !is(T == ubyte[][])){
                buffer.length = T.sizeof*fileMeta.rowGroups[0].num_rows;
            }
            else{
                buffer.length = fileMeta.rowGroups[0].num_rows;
            }
            file = f;
        }

        static if (!isRaw)
        public T fullColumn(){
            uint num_rows = reduce!((a,b) => a + b)(map!(row_group => row_group.num_rows)
                            (fileMeta.rowGroups));

            T column;
            column.length = num_rows;

            foreach(i, elem; this){
                column[i] = elem;
            }

            return column;
        }

        private static ElementType!T parse(ref byte[] bytes){
            // variable size fields require additional parsing step
            static if(is(ElementType!T==ubyte[])){
                auto length = *(cast( int*)(bytes.ptr));
                ElementType!T value = cast(ubyte[])bytes[int.sizeof..int.sizeof+length];//.dup;
                bytes = bytes[int.sizeof+length..$];
            }
            else{
                // To parse row groups
                static if(isRaw){
                    int localOffset = 0;
                }
                else static if(bfield == bamFields.bin 
                        || bfield == bamFields.flag){
                    int localOffset = 2;  
                }
                else static if(bfield == bamFields.mapq){
                    int localOffset = 1;
                }
                else{
                    int localOffset = 0;
                }
                ElementType!T value = *(cast( ElementType!T*)(bytes.ptr + localOffset));
                bytes = bytes[int.sizeof..$]; // fields smaller than int packed into one 
                //writeln(localOffset);
            }
            return value;
        }
        
        static if(!isRaw)
        int opApply(int delegate(ElementType!T) dg){
            int result = 0;
            foreach(rowGroup; fileMeta.rowGroups){
                buffer.length = rowGroup.num_rows;
                fetchBuffer(rowGroup);
                foreach(ref ElementType!T elem; buffer){
                    result = dg(elem);
                    if(result) break;
                }
            }
            
            return result;
        }

        static if(!isRaw)
        int opApply(int delegate(ref size_t i, ElementType!T) dg){
            int result = 0;
            size_t i = 0;
            foreach(rowGroup; fileMeta.rowGroups){
                buffer.length = rowGroup.num_rows;
                fetchBuffer(rowGroup);
                foreach(ElementType!T elem; buffer){
                    result = dg(i, elem);
                    ++i;
                    if(result) break;
                }
            }

            return result;
        }

        private void fetchBuffer(RowGroupMeta rowGroup){
            byte[] rawBuffer;

            ulong size = rowGroup.columnsSizes[columnType];
            rawBuffer.length = size;
            
            file.seek(rowGroup.columnsOffsets[columnType]);
            file.rawRead(rawBuffer);
            auto plainBuffer = Snappy.uncompress(rawBuffer);

            static if(isRaw && !is(T==ubyte[][])){
                buffer = plainBuffer;
            }
            else{
                foreach(i, ref elem; buffer){
                    elem = parse(plainBuffer);
                }
            }
        }

        static if(!isRaw)
        T getChunk(int i) {
            fetchBuffer(fileMeta.rowGroups[i]); 
            return buffer; 
        }
        else
        ubyte[] getChunk(int i) {
            fetchBuffer(fileMeta.rowGroups[i]); 
            return buffer; 
        }
        
        private: 
            static if(isRaw){
                ubyte[] buffer;
            }
            else{
                T buffer;
            }
            const ColumnTypes columnType = bamfToCbamf[bfield];
            FileMeta fileMeta;
            File file;
    }
    
    unittest{
        auto fn1 = getcwd() ~ "/source/tests/test3.cbam";
        auto fn2 = getcwd() ~ "/source/tests/test1.bam";
        
        File file = File(fn1, "r");
        FileReader fileR = new FileReader(file);

        auto CBAMpos  = fileR.getColumn!(int[], bamFields.pos);
        auto CBAMflag = fileR.getColumn!(short[], bamFields.flag);
        auto CBAMmapq = fileR.getColumn!(ubyte[], bamFields.mapq);
        auto CBAMqual = fileR.getColumn!(ubyte[][], bamFields.raw_qual);

        auto posCol  = CBAMpos.fullColumn();
        auto flagCol = CBAMflag.fullColumn();
        auto mapqCol = CBAMmapq.fullColumn();
        auto qualCol = CBAMqual.fullColumn();

        BamReadBlobStream reader = BamReadBlobStream(fn2);
        int i = 0;
        while(!reader.empty()){
            auto temp = reader.front();

            assert(posCol[i] == temp.pos);
            assert(flagCol[i] == temp._flag);
            assert(mapqCol[i] == temp._mapq);
            assert(equal(qualCol[i], temp.raw_qual));
            reader.popFront();
            ++i;
        }
        assert(posCol.length-1 == i);
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