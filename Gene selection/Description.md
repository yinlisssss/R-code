The original TCGA input file is larger than 25MB and exceeds the limit for uploading to Github

So, the file is split using command 'split'

```
split -b 20000kb TCGALUAD.zip TCGA_input_part
```

To recover them, simply use
```
cat TCGA_input_part* > TCGALUAD.zip
```

