# What a losslessly decoded boxplot can prove

This case study computes sharp bounds on sample means compatible with the boxplots in an actual research figure. The extractor reads vector geometry. A separate exact optimizer constructs the smallest- and largest-mean samples with the same boxplot marks.

## Result

Figure 3 from the authors’ repository for Zapata and Villarroel’s 2019 *Escallonia* paper contains 39 boxplots. The program recovered the three quartiles, both whiskers and all 30 outlier marks, including multiplicities, for all 39 species. These values agree with the summaries recomputed from the 746 nonmissing flower counts in the authors’ CSV.

Even given each true sample size, the flower-count mean is ambiguous for **27 of 39** species. Among the **738** pairs with unequal source means, **46** admit an opposite strict ordering while preserving both boxplots. This is a pairwise existence count; different reversals can use different compatible readings of the same species.

| Species | Sample size | Minimum compatible mean | Mean in source data | Maximum compatible mean |
| --- | ---: | ---: | ---: | ---: |
| E. pendula | 28 | 156.14 | 169.29 | 191.86 |
| E. pulverulenta | 20 | 165.50 | 171.50 | 174.50 |
| E. paniculata | 62 | 126.84 | 148.94 | 171.92 |
| E. herrerae | 5 | 157.00 | 157.00 | 157.00 |

An explicit reversal for the first pair is obtained by replacing *E. pendula* with one 80, seven 136s, seven 180s, seven 200s and six 280s. Replace *E. pulverulenta* with four 140s, five 150s, ten 180s and one 200. Both boxes, whiskers, outlier lists and sample sizes are preserved exactly. The source means are ordered 169.29 < 171.50, while the replacement means are 191.86 > 165.50. The JSON output gives both extreme samples for every species.

The analysis concerns what Figure 3 encodes about sample means. The authors used other analyses to support their species description. The figure-specific PDF decoder was validated on this figure, and the integer optimization applies to type-7 quartiles, known sample counts and known outlier multiplicities.

## Why the bounds are sharp

A solution is a sorted sequence of integer flower counts. The visible outliers fix the prefix and suffix. Whiskers constrain the rest to an interval, fixing its endpoints wherever a whisker lies strictly beyond the corresponding quartile. The remaining constraints come from the three type-7 quantiles. Each is a linear equation involving either one rank or two adjacent ranks.

For each rank and possible current value, dynamic programming stores the minimum or maximum attainable prefix sum. An unconstrained transition uses the best predecessor no larger than the current value. For a two-rank quantile equation, the current value determines the only possible predecessor. By induction these transitions enumerate the feasible sequences through their optimal prefixes. Backtracking recovers a sample attaining each bound. Every recovered sample is checked again against the original boxplot summary. The independent self-test compares both bounds against brute force on 10,144 small cases.

## Reproduce

```sh
python -m pip install pymupdf==1.26.6
git clone https://github.com/zapataf/ms_eharrisii.git source-data
git -C source-data checkout ac80268e79b79480dd1ee432e8fe9f37dca0276e
python study.py source-data/Figure3.pdf source-data/data.csv --self-test > results.json
```

The PDF and CSV are obtained from the authors’ repository. SHA-256 digests of the inputs used here are:

| Input | SHA-256 |
| --- | --- |
| Figure3.pdf | `9aa24d9d12aa8d0262c3c5feaf441e780f37141c5cdd17fefa57544cbd1e4576` |
| data.csv | `aeb8209663eec62d8b35da7fb225d1acdd5ad68d17d4dc517cd7cdff14dcb303` |

Sources: [Zapata and Villarroel, PeerJ 7:e6328](https://doi.org/10.7717/peerj.6328), [author repository](https://github.com/zapataf/ms_eharrisii), and [ggplot2 boxplot conventions](https://ggplot2.tidyverse.org/reference/geom_boxplot.html).
