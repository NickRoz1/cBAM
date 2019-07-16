module tests.benchmark;

import bio.std.experimental.hts.bam.reader;
import bio.std.hts.bam.reader;
import reader;

import std.datetime.stopwatch;
import std.stdio;
import std.file;
import std.parallelism;

alias bamFields = FileReader.bamFields;

void main(){
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