# helix
contains a collection of useful scripts for the knights lab

## filter_fastq.py
```
mkdir R2
cp *R2*.fastq R2/
cd R2
for file in * ; do mv -v "$file" "$(echo $file | sed 's/_/./g' | sed 's/.R2//g')"; done

python filter_fastq_test.py -f R1/<input_fastq> -b <blast_file> -o <output>

# Do stuff to files
for file in * ; do mv -v "$file" "$(echo $file | sed 's/.001.fq/.R2.001.fq/g')"; done
```
