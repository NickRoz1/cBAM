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

    void close(){
        file.close();
    }

    enum colTypeToOff = [
        "_bin_mq_nl"         : "Offset.bin_mq_nl",
        "_flag_nc"           : "Offset.flag_nc",
        "sequence_length"    : "Offset.l_seq",
        "_next_refID"        : "Offset.next_refID",
        "_next_pos"          : "Offset.next_pos",
        "_tlen"              : "Offset.tlen",
        "read_name"          : "_read_name_offset",
        "raw_cigar"          : "_cigar_offset",
        "raw_sequence"       : "_seq_offset",
        "raw_qual"           : "_qual_offset"
    ];

    RawReadBlob[] readRowGroup(const int i){
        enforce(i >= 0 && i < fileMeta.rowGroups.length, "No such rowgroup " ~ to!string(i));

        auto sizesCol = getColumn!(bamFields.blob_size).getChunk(i);
        
        RawReadBlob[] readsBuf;
        readsBuf.length = sizesCol.length;

        foreach(blob_n, ref elem; readsBuf){
            elem._data.length = sizesCol[blob_n];
        }

        template GenFiller(ColumnTypes colType){
            static if(colType == ColumnTypes._refID || colType == ColumnTypes._pos){
                enum ternar = [ColumnTypes._refID : "refid",
                               ColumnTypes._pos   : "pos"];
                const char[] GenFiller = "
                    foreach(read_n, ref blob; readsBuf){
                        blob." ~ ternar[colType] ~ " = *(cast( int*)(col.ptr + int.sizeof*read_n));
                    }";
            }
            else static if(colType == ColumnTypes.read_name ||
                           colType == ColumnTypes.raw_cigar ||
                           colType == ColumnTypes.raw_sequence ||
                           colType == ColumnTypes.raw_qual){
                enum offset = "blob." ~ colTypeToOff[to!string(colType)];               
                const char[] GenFiller = "
                    foreach(read_n, ref blob; readsBuf){
                        blob._data["~ offset ~ ".." ~ offset ~ "+col[read_n].length] = col[read_n];
                    }";              
            }
            else{
                enum offset = colTypeToOff[to!string(colType)];
                const char[] GenFiller = "
                    foreach(read_n, ref blob; readsBuf){
                        blob._data[" ~ offset ~ ".." ~ offset ~ "+int.sizeof] = col[read_n*int.sizeof..(read_n+1)*int.sizeof];
                    }";
                }
        }

        void localReader(ColumnTypes colType)(){
            auto col = getBaseColumn!(colType).getRawChunk(i);
            mixin(GenFiller!(colType));
        }
        static foreach(member; EnumMembers!ColumnTypes){
            static if(member != ColumnTypes._blob_size){
                localReader!(member);
            }
        }

        return readsBuf;
    }

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

    auto getBaseColumn(ColumnTypes colType)(){
        static if(colType == ColumnTypes._refID ||
                  colType == ColumnTypes._pos   ||
                  colType == ColumnTypes._blob_size ||
                  colType == ColumnTypes._bin_mq_nl ||
                  colType == ColumnTypes._flag_nc ||
                  colType == ColumnTypes.sequence_length ||
                  colType == ColumnTypes._next_refID ||
                  colType == ColumnTypes._next_pos ||
                  colType == ColumnTypes._tlen){
                    return new Column!(int[], bamFields.refID, colType)(fileMeta, file);
        }
        else{
            return new Column!(ubyte[][], bamFields.refID, colType)(fileMeta, file);
        }
    }

    auto getColumn(bamFields bfield)(){
        static if(bfield == bamFields.refID ||
                  bfield == bamFields.pos   ||
                  bfield == bamFields.nextRefID ||
                  bfield == bamFields.nextPos ||
                  bfield == bamFields.tlen ||
                  bfield == bamFields.blob_size){
                    return new Column!(int[], bfield)(fileMeta, file);
                }
        static if(bfield == bamFields.bin  ||
                  bfield == bamFields.flag){
                    return new Column!(short[], bfield)(fileMeta, file);
                }
        static if(bfield == bamFields.mapq)
                    return new Column!(ubyte[], bfield)(fileMeta, file);
        static if(bfield == bamFields.read_name ||
                  bfield == bamFields.cigar     ||
                  bfield == bamFields.raw_qual  ||
                  bfield == bamFields.raw_sequence){
                    return new Column!(ubyte[][], bfield)(fileMeta, file);
                }
        assert(0);
    }

    struct Column(T, bamFields bfield, ColumnTypes columnType = bamfToCbamf[bfield])
        if(is(T == ubyte[][]) || is(T == int[]) 
        || is(T == short[])   || is(T == ubyte[])) {
            
        this(FileMeta meta, File f){
            fileMeta = meta;
            buffer.length = fileMeta.rowGroups[0].num_rows;
            file = f;
        }

        public @property long size(){
            return reduce!((a,b) => a + b)(map!(row_group => row_group.num_rows)
                            (fileMeta.rowGroups));
        }

        
        public T fullColumn(){
            T column;
            column.length = this.size;

            foreach(i, elem; this){
                column[i] = elem;
            }

            return column;
        }

        private static ElementType!T parse(ref byte[] bytes, ref long offset){
            // variable size fields require additional parsing step
            static if(is(ElementType!T==ubyte[])){
                auto length = *(cast( int*)(bytes.ptr + offset));
                offset += int.sizeof;
                ElementType!T value = cast(ubyte[])bytes[offset..offset+length];
                offset += length;
            }
            else{
                static if(bfield == bamFields.bin 
                        || bfield == bamFields.flag){
                    int localOffset = 2;  
                }
                else static if(bfield == bamFields.mapq){
                    int localOffset = 1;
                }
                else{
                    int localOffset = 0;
                }
                ElementType!T value = *(cast( ElementType!T*)(bytes.ptr + offset + localOffset));
                offset += int.sizeof; // fields smaller than int packed into one 
            }
            return value;
        }
        
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

        private void fetchBuffer(bool isRaw = false)(RowGroupMeta rowGroup){
            byte[] rawInputBuffer;

            ulong size = rowGroup.columnsSizes[columnType];
            rawInputBuffer.length = size;

            file.seek(rowGroup.columnsOffsets[columnType]);
            file.rawRead(rawInputBuffer);
            auto plainBuffer = Snappy.uncompress(rawInputBuffer);

            static if(isRaw && !is(T==ubyte[][])){
                rawBuffer = cast(ubyte[])plainBuffer;
            }
            else{
                long offset = 0;
                foreach(i, ref elem; buffer){
                    elem = parse(plainBuffer, offset);
                }
            }
        }


        T getChunk(int i) {
            fetchBuffer(fileMeta.rowGroups[i]); 
            return buffer; 
        }

        auto getRawChunk(int i) {
            fetchBuffer!(true)(fileMeta.rowGroups[i]); 
            static if(columnType == ColumnTypes.read_name ||
                      columnType == ColumnTypes.raw_cigar     ||
                      columnType == ColumnTypes.raw_qual  ||
                      columnType == ColumnTypes.raw_sequence){
                return buffer;
            }
            else{
                return rawBuffer; 
            }
        }
        
        private: 
            ubyte[] rawBuffer;
            T buffer;
            FileMeta fileMeta;
            File file;
    }
    
    unittest{
        auto fn1 = getcwd() ~ "/source/tests/test3.cbam";
        auto fn2 = getcwd() ~ "/source/tests/test1.bam";
        
        File file = File(fn1, "r");
        FileReader fileR = new FileReader(file);
 
        auto CBAMpos  = fileR.getColumn!(bamFields.pos);//fileR.getColumn!(int[], bamFields.pos);
        auto CBAMflag = fileR.getColumn!(bamFields.flag);
        auto CBAMmapq = fileR.getColumn!(bamFields.mapq);
        auto CBAMname = fileR.getColumn!(bamFields.read_name);
        auto CBAMcig  = fileR.getColumn!(bamFields.cigar);
        auto CBAMseq  = fileR.getColumn!(bamFields.raw_sequence);
        auto CBAMqual = fileR.getColumn!(bamFields.raw_qual);

        auto posCol  = CBAMpos.fullColumn();
        auto flagCol = CBAMflag.fullColumn();
        auto mapqCol = CBAMmapq.fullColumn();
        auto readNam = CBAMname.fullColumn();
        auto cigCol = CBAMcig.fullColumn();
        auto seqCol = CBAMseq.fullColumn();
        auto qualCol = CBAMqual.fullColumn();

        BamReadBlobStream reader = BamReadBlobStream(fn2);
        int i = 0;
        while(!reader.empty()){
            auto temp = reader.front();

            assert(posCol[i] == temp.pos);
            assert(flagCol[i] == temp._flag);
            assert(mapqCol[i] == temp._mapq);
            assert(equal(readNam[i], temp.read_name));
            assert(equal(cigCol[i], temp.raw_cigar));
            assert(equal(seqCol[i], temp.raw_sequence));
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