# hmp3 - Helix MP3 Encoder, with some additions

This repository contains a version of the Helix MP3 Encoder with some maintenance work applied to build on modern systems. This targets mostly Linux-based machines (tested on x86-64 and ARM), but project files for Visual Studio 2015 are also included.

## On the license of the Helix MP3 Encoder

This repository contains a source code copy of the Helix MP3 Encoder in the directory [./hmp3](./hmp3). Please study [./hmp3/LICENSE.txt](./hmp3/LICENSE.txt) for licensing information.

## On the license of other files in this repository

Files in this repository outside of the [./hmp3](./hmp3) directory are licensed under the terms of the MIT License, unless noted otherwise.

## Compiling

A Makefile tested on Linux is provided in this directory. With the usual build systems installed, a simple `make` should create a `hmp3` program binary in `builds/release`, along with compiled object files.

A binary with debug symbols can be generated with `make debug`.

## Documentation

The encoder is documented over at the [HydrogenAudio Wiki page of the Helix MP3 encoder](https://wiki.hydrogenaud.io/index.php?title=Helix_MP3_Encoder).

## Use examples

* Create a ~128 kbps VBR MP3 file:
  
  `hmp3 input.wav output.mp3`

* Create a 128 kbps CBR MP3 file:
  
  `hmp3 input.wav output.mp3 -B64`

  (Note that `-B` denotes the bitrate **per channel**, thus stereo input files are being encoded with 128 kbps.)

* Create a ~185 kbps VBR MP3 file, encode frequencies above 16 kHz, with a highpass filter of 19 kHz applied:
  
  `hmp3 input.wav output.mp3 -V100 -HF2 -F19000`

## Lossless re-coding (A1)

For MPEG-1 VBR encodes (44.1, 48 and 32 kHz, stereo, joint stereo or mono), every
frame is re-coded after the encoder has made all its decisions: the cheapest
Huffman coding of the same quantized values (big_values boundary, region split,
tables, count1 table) and the cheapest scalefactor side info for the same step
sizes (scalefac_scale, preflag, global_gain, free values, scalefac_compress, and
scfsi chosen for both granules together). The decoded PCM does not change (bit for
bit with ffmpeg and mpg123; see `HMP3_A1=1` for minimp3); files get about 2%
smaller. Frames are re-coded into a second bit reservoir that never runs
emptier than the stock one, so no frame is ever larger than before. Encoding takes
about twice as long. See `hmp3/src/pub/a1pack.h`.

* `HMP3_A1=0` turns it off (output identical to the stock encoder).
* `HMP3_A1=1` also keeps every global_gain as coded, so even decoders that compute
  the step size from two floating-point factors (minimp3) output identical PCM;
  files are about 0.1% larger than with the default.
* `HMP3_A1STATS=1` prints what the re-coding changed.

  `HMP3_A1STATS=1 hmp3 input.wav output.mp3 -V65`
