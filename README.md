# cBAM

Tools for writing/reading cBAM files.

# Usage

## Writer

```d
import source.writer;

auto bamReader = BamReadBlobStream(".../file.bam"); 
File file = File(".../file.cbam", "w");
auto fileWriter = new FileWriter(file, bamReader.header);

while(!bamReader.empty()){
        // Add record in writing buffer
        fileWriter.addRecord(bamReader.front());
        bamReader.popFront();
}
fileWriter.addRecord(bamReader.front());
fileWriter.close();
```
## Reader

```d
import std.stdio;
import source.reader;

File file = File(".../file.cbam", "r");
auto fileReader = new FileReader(file);

// Choose desired field using bamFields enum
auto bamField = FileReader.bamFields.flag;
// Get corresponding column
auto flagColumn = FileReader.getColumn(bamField);
auto dupFlag = 0x400;
bool noDups = true;

// Iterate through column
foreach(flag; flagColumn){
    if(flag & dupFlag){
        noDups = false;
        break;
    }
}

short[] flagStorage;
flagStorage.length = flagColumn.size;

// Iterate through column with index
foreach(i, flag; flagColumn){
    flagStorage[i++] = flag;
}

// Get whole column
auto checkStorage = flagColumn.fullColumn();

assert(equal(flagStorage, checkStorage));

// Choose desired row group number
auto groupNumber = 1;
RawReadBlob[] rowGroup;

// Get rowgroup
rowGroup = fileReader.readRowGroup(groupNumber);

fileReader.close();
```






