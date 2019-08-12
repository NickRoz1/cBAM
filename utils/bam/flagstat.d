
/*
    This file is part of Sambamba.
    Copyright (C) 2012-2013    Artem Tarasov <lomereiter@gmail.com>
    Modified by Nick Rozinsky
    
    Sambamba is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    Sambamba is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

*/
module sambamba.flagstat;

/// port of samtools flagstat tool
import bio.std.hts.bam.reader;
import bio.std.experimental.hts.bam.reader;
import std.stdio, std.conv, std.parallelism, std.getopt, std.string;

import sambamba.utils.common.progressbar;

ulong[2] reads, pair_all, pair_good, first, second, single, pair_map, mapped,
         dup, diff_chr, diff_high, secondary, supplementary;

void flagStatCBAM(T)(T column){
    foreach (read; column) {
        mixin ReadFlags!(read);

        size_t failed = is_qc_fail ? 1 : 0;
        ++reads[failed];
        if (!is_unmapped_raw) ++mapped[failed];
        if (is_duplicate) ++dup[failed];
        if (is_secondary_alignment) {
            ++secondary[failed];
        } else if (is_supplementary) {
            ++supplementary[failed];
        } else if (is_paired) {
            ++pair_all[failed];
            if (is_proper_pair && !is_unmapped_raw) ++pair_good[failed];
            if (is_first_of_pair) ++first[failed];
            if (is_second_of_pair) ++second[failed];
            if (mate_is_unmapped && !is_unmapped_raw) ++single[failed];
            // if (!is_unmapped_raw && !mate_is_unmapped) {
            //     ++pair_map[failed];
            //     if (read.ref_id != read.mate_ref_id) {
            //         ++diff_chr[failed];
            //         if (read.mapping_quality >= 5)
            //             ++diff_high[failed];
            //     }
            // }
        }
    }
}
void computeFlagStatistics(R)(R alignments) {
    foreach (read; alignments) {
        size_t failed = read.failed_quality_control ? 1 : 0;
        ++reads[failed];
        if (!read.is_unmapped) ++mapped[failed];
        if (read.is_duplicate) ++dup[failed];
        if (read.is_secondary_alignment) {
            ++secondary[failed];
        } else if (read.is_supplementary) {
            ++supplementary[failed];
        } else if (read.is_paired) {
            ++pair_all[failed];
            if (read.proper_pair && !read.is_unmapped) ++pair_good[failed];
            if (read.is_first_of_pair) ++first[failed];
            if (read.is_second_of_pair) ++second[failed];
            if (read.mate_is_unmapped && !read.is_unmapped) ++single[failed];
            if (!read.is_unmapped && !read.mate_is_unmapped) {
                ++pair_map[failed];
                if (read.ref_id != read.mate_ref_id) {
                    ++diff_chr[failed];
                    if (read.mapping_quality >= 5)
                        ++diff_high[failed];
                }
            }
        }
    }
}

void writeParam(string description, ulong[2] param) {
    writefln("%s + %s %s", param[0], param[1], description);
}

float percent(ulong a, ulong b) { return to!float(a) / b * 100.0; }

string percentStr(ulong a, ulong b) {
    if (b == 0) return "N/A";
    return format("%.2f%%", percent(a, b));
}

void writeParamWithPercentage(string description, ulong[2] param, ulong[2] total) {
    writefln("%s + %s %s (%s:%s)", param[0], param[1], description,
             percentStr(param[0], total[0]),
             percentStr(param[1], total[1]));
}

version(standalone) {
    int main(string[] args) {
        return flagstat_main(args);
    }
}

void printUsage() {
    stderr.writeln("Usage: sambamba-flagstat [options] <input.bam>");
    stderr.writeln();
    stderr.writeln("OPTIONS: -t, --nthreads=NTHREADS");
    stderr.writeln("            use NTHREADS for decompression");
    stderr.writeln("         -p, --show-progress");
    stderr.writeln("            show progressbar in STDERR");
    stderr.writeln("         -c, --cbam");
    stderr.writeln("            input CBAM");
}

int flagstat_main(string[] args) {
    size_t threads = totalCPUs;
    bool show_progress;
    bool cbam;

    try {
        getopt(args,
               std.getopt.config.caseSensitive,
               "nthreads|t",      &threads,
               "show-progress|p", &show_progress,
               "cbam|c", &cbam);

        if (args.length < 2) {
            printUsage();
            return 1;
        }

        if(cbam){
            File file = File(".../file.cbam", "r");
            auto fileReader = new FileReader(file);
            auto bamField = FileReader.bamFields.flag;
            auto flagColumn = FileReader.getColumn(bamField);

            flagStatCBAM(flagColumn);
        }
        else{
            auto task_pool = new TaskPool(threads);
            scope(exit) task_pool.finish();

            auto bam = new BamReader(args[1], task_pool);
            bam.assumeSequentialProcessing();

            if (show_progress) {
                auto bar = new shared(ProgressBar)();
                computeFlagStatistics(bam.readsWithProgress((lazy float p) { bar.update(p); }));
                bar.finish();
            } else {
                computeFlagStatistics(bam.reads!withoutOffsets);
            }
        }

        scope(exit) {
            writeParam("in total (QC-passed reads + QC-failed reads)", reads);
            writeParam("secondary", secondary);
            writeParam("supplementary", supplementary);
            writeParam("duplicates", dup);
            writeParamWithPercentage("mapped", mapped, reads);
            writeParam("paired in sequencing", pair_all);
            writeParam("read1", first);
            writeParam("read2", second);
            writeParamWithPercentage("properly paired", pair_good, pair_all);
            writeParam("with itself and mate mapped", pair_map);
            writeParamWithPercentage("singletons", single, pair_all);
            writeParam("with mate mapped to a different chr", diff_chr);
            writeParam("with mate mapped to a different chr (mapQ>=5)", diff_high);
        }
    } catch (Throwable e) {
        stderr.writeln(e.msg);
    }
    return 0;
}
