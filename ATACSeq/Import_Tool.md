`universalmotif` is a powerful tool to import most common motif types into R, and convert the motif information between different classes as defined by other Bioconductor packages. 
Reference: [universalmotif](https://github.com/bjmt/universalmotif)

```
library(universalmotif)
meme_motif <- read_meme("./data/mouse.meme")
meme_motif
meme_motif_pwm <- convert_motifs(meme_motif, "TFBSTools-PWMatrix")
meme_motif_pwm
```
