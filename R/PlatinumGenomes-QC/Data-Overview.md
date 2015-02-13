<!-- R Markdown Documentation, DO NOT EDIT THE PLAIN MARKDOWN VERSION OF THIS FILE -->

<!-- Copyright 2015 Google Inc. All rights reserved. -->

<!-- Licensed under the Apache License, Version 2.0 (the "License"); -->
<!-- you may not use this file except in compliance with the License. -->
<!-- You may obtain a copy of the License at -->

<!--     http://www.apache.org/licenses/LICENSE-2.0 -->

<!-- Unless required by applicable law or agreed to in writing, software -->
<!-- distributed under the License is distributed on an "AS IS" BASIS, -->
<!-- WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. -->
<!-- See the License for the specific language governing permissions and -->
<!-- limitations under the License. -->

# Part 1: Data Overview

The following example makes use of [Illumina Platinum Genomes](http://www.illumina.com/platinumgenomes/).  For more detail about how this data was loaded into the Google Genomics API, please see [Google Genomics Public Data](https://cloud.google.com/genomics/data/platinum-genomes).






```r
# By default this codelab runs upon the Illumina Platinum Genomes Variants.
# Change the table here if you wish to run these queries against a different table.
#tableReplacement <- list("_THE_TABLE_"="genomics-public-data:platinum_genomes.variants")

tableReplacement <- list("_THE_TABLE_"="google.com:biggene:test.pgp_variants_20150205")
```

## Variants

Let's take a look at a few of the [variants within BRCA1 via BigQuery](https://github.com/googlegenomics/getting-started-bigquery/blob/master/RMarkdown/literate-programming-demo.md#data-visualization):

```r
result <- DisplayAndDispatchQuery("https://raw.githubusercontent.com/googlegenomics/getting-started-bigquery/master/sql/variant-level-data-for-brca1.sql",
                                  project=project,
                                  replacements=tableReplacement)
```

```
# Retrieve variant-level information for BRCA1 variants.
SELECT
  reference_name,
  start,
  end,
  reference_bases,
  GROUP_CONCAT(alternate_bases) WITHIN RECORD AS alternate_bases,
  quality,
  GROUP_CONCAT(filter) WITHIN RECORD AS filter,
  GROUP_CONCAT(names) WITHIN RECORD AS names,
  COUNT(call.call_set_name) WITHIN RECORD AS num_samples,
FROM
  [google.com:biggene:test.pgp_variants_20150205]
WHERE
  reference_name = 'chr17'
  AND start BETWEEN 41196311
  AND 41277499
OMIT RECORD IF EVERY(alternate_bases IS NULL)
ORDER BY
  start,
  alternate_bases
Running query:   RUNNING  2.6s
```
Number of rows returned by this query: 1186.

Displaying the first few rows of the dataframe of results:
<!-- html table generated in R 3.1.1 by xtable 1.7-4 package -->
<!-- Fri Feb 13 10:30:38 2015 -->
<table border=1>
<tr> <th> reference_name </th> <th> start </th> <th> end </th> <th> reference_bases </th> <th> alternate_bases </th> <th> quality </th> <th> filter </th> <th> names </th> <th> num_samples </th>  </tr>
  <tr> <td> chr17 </td> <td align="right"> 41196362 </td> <td align="right"> 41196363 </td> <td> C </td> <td> T </td> <td align="right"> 0.00 </td> <td>  </td> <td> dbsnp.117:rs8176320,dbsnp.117:rs8176320 </td> <td align="right">   6 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41196407 </td> <td align="right"> 41196408 </td> <td> G </td> <td> A </td> <td align="right"> 0.00 </td> <td>  </td> <td> dbsnp.52:rs12516,dbsnp.52:rs12516 </td> <td align="right"> 112 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41196407 </td> <td align="right"> 41196407 </td> <td>  </td> <td> G </td> <td align="right"> 0.00 </td> <td>  </td> <td>  </td> <td align="right">   1 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41196821 </td> <td align="right"> 41196823 </td> <td> TT </td> <td> ?T </td> <td align="right"> 0.00 </td> <td>  </td> <td> dbsnp.126:rs33947868;dbsnp.129:rs60038333;dbsnp.130:rs68017638,dbsnp.126:rs33947868;dbsnp.129:rs60038333;dbsnp.130:rs68017638 </td> <td align="right">   4 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41196821 </td> <td align="right"> 41196823 </td> <td> TT </td> <td> ?TT </td> <td align="right"> 0.00 </td> <td>  </td> <td> dbsnp.126:rs33947868;dbsnp.129:rs60038333;dbsnp.130:rs68017638 </td> <td align="right">   3 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41196821 </td> <td align="right"> 41196823 </td> <td> TT </td> <td> T </td> <td align="right"> 0.00 </td> <td>  </td> <td> dbsnp.126:rs33947868;dbsnp.129:rs60038333;dbsnp.130:rs68017638,dbsnp.126:rs33947868;dbsnp.129:rs59541324;dbsnp.129:rs60038333;dbsnp.130:rs68017638,dbsnp.126:rs33947868;dbsnp.129:rs60038333;dbsnp.130:rs68017638,dbsnp.126:rs33947868;dbsnp.129:rs59541324;dbsnp.129:rs60038333;dbsnp.130:rs68017638 </td> <td align="right">   6 </td> </tr>
   </table>

## Non-Variant Segments
The source Platinum Genomes data loaded into the Google Genomics API was in [genome VCF](https://sites.google.com/site/gvcftools/home/about-gvcf/gvcf-conventions) (gVCF) format.

Let's take a look at a few non-variant segments within BRCA1:

```r
result <- DisplayAndDispatchQuery("./sql/non-variant-segments.sql",
                                  project=project,
                                  replacements=tableReplacement)
```

```
# Retrieve non-variant segments for BRCA1, flattening by sample.
SELECT
  reference_name,
  start,
  end,
  reference_bases,
  GROUP_CONCAT(alternate_bases) WITHIN RECORD AS alternate_bases,
  call.call_set_name,
  GROUP_CONCAT(STRING(call.genotype)) WITHIN call AS genotype,
FROM
  [google.com:biggene:test.pgp_variants_20150205]
WHERE
  reference_name = 'chr17'
  AND start BETWEEN 41196311
  AND 41277499
OMIT RECORD IF SOME(alternate_bases IS NOT NULL)
ORDER BY
  start,
  call.call_set_name
Running query:   RUNNING  2.3sRunning query:   RUNNING  3.0sRunning query:   RUNNING  3.6s
Retrieving data:  7.3sRetrieving data: 11.7sRetrieving data: 15.8sRetrieving data: 20.8sRetrieving data: 26.6sRetrieving data: 30.6s
```
Number of rows returned by this query: 63694.

Displaying the first few rows of the dataframe of results:
<!-- html table generated in R 3.1.1 by xtable 1.7-4 package -->
<!-- Fri Feb 13 10:31:15 2015 -->
<table border=1>
<tr> <th> reference_name </th> <th> start </th> <th> end </th> <th> reference_bases </th> <th> alternate_bases </th> <th> call_call_set_name </th> <th> genotype </th>  </tr>
  <tr> <td> chr17 </td> <td align="right"> 41196363 </td> <td align="right"> 41196821 </td> <td> = </td> <td>  </td> <td> hu0E64A1 </td> <td> 0,0 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41196363 </td> <td align="right"> 41196821 </td> <td> = </td> <td>  </td> <td> hu3A8D13 </td> <td> 0,0 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41196363 </td> <td align="right"> 41196817 </td> <td> = </td> <td>  </td> <td> hu553620 </td> <td> 0,0 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41196363 </td> <td align="right"> 41196407 </td> <td> = </td> <td>  </td> <td> huA4F281 </td> <td> 0,0 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41196363 </td> <td align="right"> 41196821 </td> <td> = </td> <td>  </td> <td> huEBD467 </td> <td> 0,0 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41196363 </td> <td align="right"> 41196821 </td> <td> = </td> <td>  </td> <td> huFCC1C1 </td> <td> 0,0 </td> </tr>
   </table>
So for any analyses that require us to know for example _"how many samples do and do not have a particular SNP?"_, we'll need to make sure that the non-variant segments are considered in addition to the variants.

Note that Complete Genomics data also includes non-variant segments and requires the same consideration.

## Alternative Allele Field

And then let's take a look at the domain and range of values for alternate_bases:

```r
result <- DisplayAndDispatchQuery("./sql/characterize-alts.sql",
                                  project=project,
                                  replacements=tableReplacement)
```

```
# Variants are only SNPs and INDELs, with no special characters.
SELECT
  COUNT(1) AS number_of_variant_records,
  REGEXP_MATCH(alternate_bases,
    r'^[A,C,G,T]+$') AS alt_contains_no_special_characters,
  MAX(LENGTH(reference_bases)) AS max_ref_len,
  MAX(LENGTH(alternate_bases)) AS max_alt_len
FROM
  [google.com:biggene:test.pgp_variants_20150205]
OMIT RECORD IF EVERY(alternate_bases IS NULL)
GROUP BY
  alt_contains_no_special_characters
```

<!-- html table generated in R 3.1.1 by xtable 1.7-4 package -->
<!-- Fri Feb 13 10:31:18 2015 -->
<table border=1>
<tr> <th> number_of_variant_records </th> <th> alt_contains_no_special_characters </th> <th> max_ref_len </th> <th> max_alt_len </th>  </tr>
  <tr> <td align="right"> 40284485 </td> <td> TRUE </td> <td align="right"> 198 </td> <td align="right"> 215 </td> </tr>
  <tr> <td align="right"> 2070037 </td> <td> FALSE </td> <td align="right"> 200 </td> <td align="right"> 191 </td> </tr>
   </table>
We see from the query results that there are no special charaters in alternate_bases and the maximum length is ~50 base pairs.

## Genotype Field

And finally let's take a look at the domain and range of values for genotype:

```r
result <- DisplayAndDispatchQuery("./sql/genotypes-brca1.sql",
                                  project=project,
                                  replacements=tableReplacement)
```

```
# Query to show the variety of genotypes within BRCA1, such as
# single allele genotypes.
SELECT
  genotype,
  COUNT(genotype) AS genotype_count
FROM (
  SELECT
    GROUP_CONCAT(STRING(call.genotype)) WITHIN call AS genotype,
  FROM
  [google.com:biggene:test.pgp_variants_20150205]
  WHERE
    reference_name = 'chr17'
    AND start BETWEEN 41196311
    AND 41277499
    )
GROUP BY
  genotype
ORDER BY
  genotype_count DESC
```

<!-- html table generated in R 3.1.1 by xtable 1.7-4 package -->
<!-- Fri Feb 13 10:31:20 2015 -->
<table border=1>
<tr> <th> genotype </th> <th> genotype_count </th>  </tr>
  <tr> <td> 0,0 </td> <td align="right"> 40485 </td> </tr>
  <tr> <td> -1,-1 </td> <td align="right"> 21781 </td> </tr>
  <tr> <td> 1,0 </td> <td align="right"> 10908 </td> </tr>
  <tr> <td> 1,1 </td> <td align="right"> 4152 </td> </tr>
  <tr> <td> 1,-1 </td> <td align="right"> 2093 </td> </tr>
  <tr> <td> -1,0 </td> <td align="right"> 1008 </td> </tr>
  <tr> <td> 0,-1 </td> <td align="right"> 420 </td> </tr>
  <tr> <td> -1,1 </td> <td align="right">  87 </td> </tr>
  <tr> <td> 1,2 </td> <td align="right">  60 </td> </tr>
  <tr> <td> 0,1 </td> <td align="right">   7 </td> </tr>
   </table>
We see from the query results the variety of genotypes within BRCA1.

# Summary

To summarize attributes of this particular dataset that we need to consider when working with Platinum Genomes data:
* It has non-variant segments which adds complexity above and beyond [similar examples for the 1,000 Genomes dataset](https://github.com/googlegenomics/bigquery-examples/blob/master/1000genomes/sql/README.md).
* It is comprised only of SNPs and INDELs (contains no structural variants).
* The values for `alternate_bases` are just comprised of the letters A,C,G,T (e.g., contains no `<DEL>` values).
* It contains some single-allele and 1/2 genotypes.

--------------------------------------------------------
_Next_: [Part 2: Data Conversion](./Data-Conversion.md)
