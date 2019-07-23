module tests.benchmark;

//import bio.std.experimental.hts.bam.reader;
import bio.std.hts.bam.reader;
import reader;
import utils.bam.reader;

import std.datetime.stopwatch;
import std.stdio;
import std.file;
import std.parallelism;

alias bamFields = FileReader.bamFields;

void main1(){
    //posColumn();
    qualColumn();
}

void posColumn(){
    auto fn1 = getcwd() ~ "/source/tests/test5.cbam";
    File file = File(fn1, "r");
    FileReader cbam_reader = new FileReader(file);

    auto fn2 = getcwd() ~ "/source/tests/test2.bam";
    BamReadBlobs bam_reader_exp = BamReadBlobs(fn2);

    auto sw = StopWatch(AutoStart.no);
    sw.start();
    int[] cbam_test = cbam_reader.getColumn!(int[])(bamFields.pos);
    sw.stop();
    auto cbam_dur = sw.peek.total!"usecs";
    sw.reset();

    int[] bam_exp_test;
    bam_exp_test.length = cbam_test.length;
    int idx = 0;
    sw.start();
    foreach(blob; bam_reader_exp){
        bam_exp_test[idx] = blob.pos;
        ++idx;
    }
    sw.stop();
    auto bam_exp_dur = sw.peek.total!"usecs";
    sw.reset();

    
    auto pool = new TaskPool(1);
    scope (exit) pool.finish();
    BamReader bam_reader = new BamReader(fn2, pool);
    int[] bam_test;
    bam_test.length = cbam_test.length;
    idx = 0;
    sw.start();
    foreach(read; bam_reader.reads){
        bam_test[idx] = read.position;
        ++idx;
    }
    sw.stop();
    auto bam_dur = sw.peek.total!"usecs";
    sw.reset();

    writeln("CBAM reader pos column parsing duration: ", cbam_dur);
    writeln("Experimental BAM reader pos column parsing duration: ", bam_exp_dur);
    writeln("BAM reader pos column parsing duration: ", bam_dur);
    writefln("CBAM reader %s times faster than experimental BAM reader", cast(double)bam_exp_dur/cbam_dur);
    writefln("CBAM reader %s times faster than BAM reader", cast(double)bam_dur/cbam_dur); 
}

void qualColumn(){
    auto fn1 = getcwd() ~ "/source/tests/test5.cbam";
    File file = File(fn1, "r");
    FileReader cbam_reader = new FileReader(file);

    auto fn2 = getcwd() ~ "/source/tests/test2.bam";
    BamReadBlobs bam_reader_exp = BamReadBlobs(fn2);

    auto sw = StopWatch(AutoStart.no);
    sw.start();
    ubyte[][] cbam_test = cbam_reader.getColumn!(ubyte[][])(bamFields.raw_qual);
    sw.stop();
    auto cbam_dur = sw.peek.total!"usecs";
    sw.reset();

    ubyte[][] bam_exp_test;
    bam_exp_test.length = cbam_test.length;
    int idx = 0;
    sw.start();
    foreach(blob; bam_reader_exp){
        bam_exp_test[idx] = blob.raw_qual;
        ++idx;
    }
    sw.stop();
    auto bam_exp_dur = sw.peek.total!"usecs";
    sw.reset();

    
    auto pool = new TaskPool(1);
    scope (exit) pool.finish();
    BamReader bam_reader = new BamReader(fn2, pool);
    ubyte[][] bam_test;
    bam_test.length = cbam_test.length;
    idx = 0;
    sw.start();
    foreach(read; bam_reader.reads){
        bam_test[idx] = read.base_qualities;
        ++idx;
    }
    sw.stop();
    auto bam_dur = sw.peek.total!"usecs";
    sw.reset();

    writeln("CBAM reader raw_qual column parsing duration: ", cbam_dur);
    writeln("Experimental BAM reader raw_qual column parsing duration: ", bam_exp_dur);
    writeln("BAM reader raw_qual column parsing duration: ", bam_dur);
    writefln("CBAM reader %s times faster than experimental BAM reader", cast(double)bam_exp_dur/cbam_dur);
    writefln("CBAM reader %s times faster than BAM reader", cast(double)bam_dur/cbam_dur); 
}

const string dir = "/";
void bench(string args[]){
    if(args.length == 0){
        writeln("No args");
        return;
    }
    if(args[1] == "convert"){
        if(args.length != 4){
            writeln("Wrong input");
            return;
        }
        auto fn = getcwd() ~ dir ~ args[2];
        auto ofn = getcwd() ~ dir ~ args[3];
        writeln(bamToCbam(fn, ofn));
        return;
    }
    if(args[1] == "bench"){
        if(args.length != 3){
            writeln("Wrong input {kek}");
            return;
        }
        auto fn = getcwd() ~ dir ~ args[2];
        
        
        version(CBAM){
            File file = File(fn, "r");
            FileReader cbam_reader = new FileReader(file);

            void get_fixed_field(){
                auto perf = cbam_reader.getColumn!(int[])(FileReader.bamFields.pos);
            }
            void get_variable_field(){
                 auto perf = cbam_reader.getColumn!(ubyte[][])(FileReader.bamFields.raw_qual);
            }
            void whole_file(){
                foreach(int i; 0..cast(int)cbam_reader.fileMeta.rowGroups.length){
                    cbam_reader.readRowGroup(i);
                }
            }
        }
        version(BioDexp){           
            
            void get_fixed_field(){
                int[] buffer;
                auto exp_reader = BamReadBlobs(fn);

                foreach(read; exp_reader){
                    buffer ~= read.refid;
                }
            }
            void get_variable_field(){
                ubyte[][] buffer;
                auto exp_reader = BamReadBlobs(fn);

                foreach(read; exp_reader){
                    buffer ~= read.raw_qual;
                }
            }
            void whole_file(){
                auto exp_reader = BamReadBlobs(fn);

                foreach(read; exp_reader){
                    //
                }
            }
        }
        version(BioDlegacy){           
            auto pool = new TaskPool(1); 
            scope (exit) pool.finish();  
            
            void get_fixed_field(){
                int[] buffer;
                auto reader = new BamReader(fn);

                foreach(read; reader.reads){
                    buffer ~= read.ref_id;
                }
            }
            void get_variable_field(){
                ubyte[][] buffer;
                auto reader = new BamReader(fn);

                foreach(read; reader.reads){
                    buffer ~= read.base_qualities;
                }
            }
            void whole_file(){
                auto reader = new BamReader(fn);

                foreach(read; reader.reads){
                    //
                }
            }
        }
        const int N = 1;
        version(CBAM){
            auto r = benchmark!(get_fixed_field,get_variable_field,whole_file)(N);
            writefln("CBAM\nFixed size field parsing took: %s\nVariable size field parsing took: %s\nWhole file parsing took: %s",
                r[0].total!"msecs"/N, r[1].total!"msecs"/N, r[2].total!"msecs"/N);
        }
        version(BioDexp){
            auto r = benchmark!(get_fixed_field,get_variable_field,whole_file)(N);
            writefln("BioD experimental reader\nFixed size field parsing took: %s\nVariable size field parsing took: %s\nWhole file parsing took: %s",
                r[0].total!"msecs"/N, r[1].total!"msecs"/N, r[2].total!"msecs"/N);
        }
        version(BioDlegacy){
            auto r = benchmark!(get_fixed_field,get_variable_field,whole_file)(N);
            writefln("BioD legacy reader\nFixed size field parsing took: %s\nVariable size field parsing took: %s\nWhole file parsing took: %s",
                r[0].total!"msecs"/N, r[1].total!"msecs"/N, r[2].total!"msecs"/N);
        }
        return;
    }
    writeln("Wrong input");
}