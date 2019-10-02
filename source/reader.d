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
import std.meta;

import source.writer;
import bio.std.experimental.hts.bam.header;
import utils.bam.reader;
import snappy.snappy;



class FileReader {
    File file;
    FileMeta fileMeta;
    BamHeader bamHeader;
    ///
    this(File f) {
        file = f;
        parse_meta();
    }

    ~this(){
        if(!file.isOpen()) return;
        file.close();
    }
    ///
    void close(){
        file.close();
    }

    /// Maps column types to Offsets
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


    /// Returns array of RawReadBlobs representing rowGroup.
    RawReadBlob[] readRowGroup(const int i){
        enforce(i >= 0 && i < fileMeta.rowGroups.length, "No such rowgroup " ~ to!string(i));

        auto sizesCol = getColumn!(bamFields.blob_size).getChunk(i);
        
        RawReadBlob[] readsBuf;
        readsBuf.length = sizesCol.length;

        /// Preallocates blob size to avoid subsequent allocations
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

    /// enum of BAM fields held in CBAM file
    enum bamFields{refID, pos, bin, mapq, flag, nextRefID,
        nextPos, tlen, read_name, cigar, raw_sequence, raw_qual, blob_size}

    /// Maps BAM fields to Columns (some BAM fields packed together)
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

    /**
        Generates Column given wanted bamField.

        Alias `colT` contains the type

        Example:
        ---
        FileReader cbam_reader = new FileReader(file);

        cbam_reader.getColType!(bfield).colT(fileMeta, file);
        ---
    */

    template getColType(bamFields bfield){
        static if(bfield == bamFields.refID ||
                  bfield == bamFields.pos   ||
                  bfield == bamFields.nextRefID ||
                  bfield == bamFields.nextPos ||
                  bfield == bamFields.tlen ||
                  bfield == bamFields.blob_size){
                    alias colT = Column!(int[], bfield);
                }
        else static if(bfield == bamFields.bin  ||
                  bfield == bamFields.flag){
                    alias colT = Column!(short[], bfield);
                }
        else static if(bfield == bamFields.mapq)
                    alias colT = Column!(ubyte[], bfield);
        else static if(bfield == bamFields.read_name ||
                  bfield == bamFields.cigar     ||
                  bfield == bamFields.raw_qual  ||
                  bfield == bamFields.raw_sequence){
                    alias colT = Column!(ubyte[][], bfield);
        }
        else static assert(false, "This BAM field is not supported yet.");
    }

    struct VariadicReader(Cols...) if(Cols.length > 0){
        
        template genSeq(Args...){
            static if (Args.length == 0){
                alias genSeq = AliasSeq!();
            }
            else{
                alias genSeq = AliasSeq!(getColType!(Args[0]).colT, genSeq!(Args[1..$]));
            }
        }
        template extractType(Args...){
            static if (Args.length == 0){
                alias extractType = AliasSeq!();
            }
            else{
                alias extractType = AliasSeq!(ElementType!(TemplateArgsOf!(Args[0])[0]), extractType!(Args[1..$]));
            }
        }

        alias colsT = genSeq!(Cols);
        colsT columns;
        alias readT = extractType!(colsT);

        this(FileMeta meta, File f){
            foreach(i, ref col; columns){
                col = colsT[i](meta, f);
            }
        }

        void getFields(ref readT fields){
            foreach(i, ref field; fields){
                field = columns[i].front;
                columns[i].popFront();
            }
        }

        int opApply(int delegate(readT) dg){
            int result = 0;

            import std.range;

            readT fields;
            for(getFields(fields); !columns[0].empty; getFields(fields)){
                result = dg(fields);
                if(result) break;
            }
        
            return result;
        }
    }
    unittest{
        auto fn1 = getcwd() ~ "/source/tests/test3.cbam";
        auto fn2 = getcwd() ~ "/source/tests/test1.bam";

        File file = File(fn1, "r");
        FileReader fileR = new FileReader(file);
        BamReadBlobStream BAMreader = BamReadBlobStream(fn2);

        auto test = VariadicReader!(bamFields.bin, bamFields.flag, bamFields.raw_qual)(fileR.fileMeta, file);
        foreach(a, b, c; test){
            auto temp = BAMreader.front();

            assert(a==temp._bin);
            assert(b==temp._flag);
            assert(c==temp.raw_qual);
            BAMreader.popFront();
        }
    }
    
    /// Get field specified by enum ColumnTypes. 
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
                      // refID as a default field, simple plug
                    return new Column!(int[], bamFields.refID, colType)(fileMeta, file);
        }
        else{ // FIX: Passes the case when user sent wrong colType
            return new Column!(ubyte[][], bamFields.refID, colType)(fileMeta, file);
        }
    }

    /**
        Returns constructed Column struct for reading requested BAM field. 

        Example:
        ---
        FileReader cbam_reader = new FileReader(file);

        auto colRefID = cbam_reader.getColumn(bamFields.refID);
        ---
    */
    auto getColumn(bamFields bfield)(){
        return getColType!(bfield).colT(fileMeta, file);
    }

    /**
        Column struct.

        Gives access to CBAM columns.

        Examples:
        ---
        
        FileReader cbam_reader = new FileReader(file);

        auto col = cbam_reader.getColumn(bamFields.refID);

        /// Print refID of every BAM record
        foreach(elem; col){
            writeln(elem);
        }

        /// With index
        foreach(i, elem; col){
            writefln("%s: %s", i, elem);
        }

        /// Get array of refID of every BAM record
        auto values = col.fullColumn();
        
        ---
    */
    struct Column(T, bamFields bfield, ColumnTypes columnType = bamfToCbamf[bfield])
        if(is(T == ubyte[][]) || is(T == int[]) 
        || is(T == short[])   || is(T == ubyte[])) {
            
        this(FileMeta meta, File f){
            enforce(f.isOpen());
            fileMeta = meta;
            buffer.length = fileMeta.rowGroups[0].num_rows;
            file = f;
            enforce(buffer.length != 0, ("File is empty"));
            rangeStatus = RangeStatus(meta.rowGroups);
            fetchBuffer(fileMeta.rowGroups[0]);
            popFront();
        }

        /// Returns number of records in whole CBAM file
        public @property long size(){
            return reduce!((a,b) => a + b)(map!(row_group => row_group.num_rows)
                            (fileMeta.rowGroups));
        }

        /// Returns array representing whole CBAM column
        public T fullColumn(){
            T column;
            column.length = this.size;

            foreach(i, elem; this){
                column[i] = elem;
            }

            return column;
        }

        /// Parses value from byte buffer
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
        
        @property ElementType!T front(){
            return rangeStatus.frontElem;
        }

        @property bool empty(){
            return rangeStatus.empty;
        }

        void popFront(){
            enforce(!rangeStatus.empty);
            
            if(rangeStatus.bufOver){
                rangeStatus.advance();
                fetchBuffer(fileMeta.rowGroups[rangeStatus.rowGroupNum]);
            }

            rangeStatus.frontElem = buffer[rangeStatus.curElem];
            rangeStatus.curElem++;
        }

        /// Iterates over CBAM column
        int opApply(int delegate(ElementType!T) dg){
            int result = 0;
            foreach(rowGroup; fileMeta.rowGroups){
                fetchBuffer(rowGroup);
                foreach(ref ElementType!T elem; buffer){
                    result = dg(elem);
                    if(result) break;
                }
            }
            
            return result;
        }

        /// Iterates over CBAM column with index
        int opApply(int delegate(ref size_t i, ElementType!T) dg){
            int result = 0;
            size_t i = 0;
            foreach(rowGroup; fileMeta.rowGroups){
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
            enforce(file.isOpen());

            buffer.length = rowGroup.num_rows;
            byte[] rawInputBuffer;

            ulong size = rowGroup.columnsSizes[columnType];
            rawInputBuffer.length = size;

            file.seek(rowGroup.columnsOffsets[columnType]);
            file.rawRead(rawInputBuffer);
            auto plainBuffer = Snappy.uncompress(rawInputBuffer);

            static if(isRaw && !is (T==ubyte[][])){
                rawBuffer = cast(ubyte[])plainBuffer;
            }
            else {

                void fillBuffer(byte[] plainBuf, T buf)
                {
                    long offset = 0;
                    foreach (i, ref elem; buf)
                    {
                        elem = parse(plainBuf, offset);
                    }
                }

                struct inputWrapper
                {
                    byte[] bytes;
                    long offset;
                }

                void fillIntegralBuffer(inputWrapper wrapper, shared(T) buf)
                {
                    import core.atomic : atomicStore;

                    long offset = 0;
                    while (offset < wrapper.bytes.length)
                    {
                        // atomicStore(buf[wrapper.offset++],
                        //         cast(shared ElementType!T) parse(wrapper.bytes, offset));
                        buf[wrapper.offset++] = cast(shared ElementType!T) parse(wrapper.bytes,
                                offset);
                    }
                }

                if (is(T == ubyte[][]))
                    fillBuffer(plainBuffer, buffer);

                else
                {
                    auto sw = StopWatch(AutoStart.no);
                    sw.start();
                    auto pool = new TaskPool();
                    scope (exit)
                        pool.finish();

                    auto chunkSize = plainBuffer.length / pool.size();
                    inputWrapper[] inputWrappers;
                    if (chunkSize * pool.size() == plainBuffer.length && chunkSize % 64 == 0)
                    {
                        inputWrappers.length = pool.size();
                    }
                    else
                    {
                        inputWrappers.length = pool.size() + 1;
                    }

                    while (chunkSize % 64 != 0)
                        chunkSize--;

                    foreach (i, ref wrapper; inputWrappers)
                    {
                        wrapper.offset = i * chunkSize / 4; // fields size is 4
                        size_t offset = i * chunkSize;

                        if (offset + chunkSize > plainBuffer.length)
                        {
                            wrapper.bytes = plainBuffer[offset .. $];
                        }
                        else
                        {
                            wrapper.bytes = plainBuffer[offset .. offset + chunkSize];
                        }
                    }
                    auto tempShared = cast(shared T) buffer;

                    foreach (wrapper; parallel(inputWrappers))
                    {
                        fillIntegralBuffer(wrapper, tempShared);
                    }
                    sw.stop();
                    //writeln("BUFFER FETCH ELAPSED: ", sw.peek.total!"msecs");
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

        string toString() {
            return format("COLUMN TYPE   : %s\nCOLUMN bfield : %s\nCOLUMN colType: %s\nCOLUMN fileno : %s\n", 
            typeid(buffer), to!string(bfield), to!string(columnType), file.fileno);
        }
        
        private: 
            ubyte[] rawBuffer;
            T buffer;

            struct RangeStatus{
                this(RowGroupMeta[] groups){ 
                    rowGroups = groups; 
                    rowGroupNum = 0;
                    curElem = 0;
                }
                RowGroupMeta[] rowGroups;
                ElementType!T frontElem;
                int rowGroupNum;
                int curElem;
                @property bool empty(){ return bufOver && exhausted; }
                @property bool bufOver(){ return curElem == rowGroups[rowGroupNum].num_rows; }
                @property bool exhausted(){ return rowGroupNum == rowGroups.length-1; }
                void advance(){
                    rowGroupNum++;
                    curElem = 0;
                }
            }
            RangeStatus rangeStatus;
            
            FileMeta fileMeta;
            public File file;
            public int status;
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

        auto rangeFlagCol = fileR.getColumn!(bamFields.flag);
        while(!reader.empty()){
            auto temp = reader.front();

            assert(rangeFlagCol.front() == temp._flag);
            if(!rangeFlagCol.empty) rangeFlagCol.popFront();

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