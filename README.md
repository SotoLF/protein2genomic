# Protein2Genomic

High-performance C++ tool for mapping protein domain coordinates to genomic coordinates using GTF annotations.

## Features

- **Ultra-fast processing**: Optimized for 100,000+ protein domains
- **Memory efficient**: Smart indexing with minimal RAM usage
- **Multi-threaded**: OpenMP support for parallel processing
- **Binary indexing**: Build once, use many times
- **Multiple output formats**: BED and detailed formats

## Building
```bash
mkdir build && cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make -j$(nproc)
