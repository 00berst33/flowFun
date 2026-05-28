# Setting up differential analysis

In this article we focus only on the steps for preparing differential
analysis, including how to specify sample information and define the
comparisons you’d like to make.

### Specifying table of sample information

Once an appropriate `GatingSet` has been prepared, performing a basic
differential analysis requires very little input from the user; however,
it is crucial that this input is correct. In order to make comparisons
between groups, R needs information about each sample in your
experiment, and what group(s) they belong to. This information can then
be used to generate a contrasts matrix and perform the tests.

To facilitate this, the user must create a table that meets a few
requirements. The first is that it contains at least 3 columns:

- `sample.name`, specifying the name of each sample.
- `filename`, specifying the file name of each sample. These names
  should correspond with the names you see when you call `sampleNames()`
  on your `GatingSet` or `flowSet`.
- A column specifying the experimental group to test within for
  differences. For example, `sex`, or `treatment`.

The second is that each sample contained in your `GatingSet` must be
represented in the table. Meaning, if you have 15 samples in your
`GatingSet`, you should have 15 corresponding rows in your table.

This table may be created programmatically in R, or exported from Excel
as a CSV file. The figure below is an example of an appropriate table
for an experiment with 4 samples.

#### Example 1

| sample.name | filename     | treatment   |
|:------------|:-------------|:------------|
| Sample 1    | sample_1.fcs | Treatment 1 |
| Sample 2    | sample_2.fcs | Treatment 2 |
| Sample 3    | sample_3.fcs | Treatment 1 |
| Sample 4    | sample_4.fcs | Treatment 2 |

This basic example would allow for testing between `"Treatment 1"` and
`"Treatment 2"`. Most experiments will have more information about their
samples, and call for more particular comparisons. As long as the
requirements listed above are adhered to, the user can add as many
columns to this table as they like. The example below shows another
example of an appropriate table, this time with more details.

#### Example 2

| sample.name | filename     | treatment | age.group | sex    | pregnancy    |
|:------------|:-------------|:----------|:----------|:-------|:-------------|
| Sample 1    | sample_1.fcs | Treatment | \<65      | female | pregnant     |
| Sample 2    | sample_2.fcs | Treatment | 65+       | female | not.pregnant |
| Sample 3    | sample_3.fcs | Ctrl      | 65+       | male   | X            |
| Sample 4    | sample_4.fcs | Treatment | 65+       | male   | X            |
| Sample 5    | sample_5.fcs | Ctrl      | \<65      | female | not.pregnant |
| Sample 6    | sample_6.fcs | Treatment | \<65      | male   | X            |

This table corresponds to an experiment with 6 samples, and has 3 more
columns than the previous example: `age.group`, `sex`, and `pregnancy`.
Note that in the column `pregnancy`, possible values are `"pregnant"`,
`"not.pregnant"`, or, when `sex == "male"`, `"X"`. Use an `"X"` (or
leave the cell empty, if making the table in Excel) when a category is
not applicable to a sample.

Many more comparisons can be created with this table, e.g. males
vs. females within the control group only, or pregnant vs. non-pregnant
females within the treatment group (of course these comparisons wouldn’t
be especially helpful with so few samples, but this table is for
demonstration purposes only).

### Define comparisons object

The other essential input expected from the user is a list defining the
comparisons to make between groups. A summary of how to define this
object can be found in the examples for the function
[`prepareSampleInfo()`](https://00berst33.github.io/flowFun/reference/prepareSampleInfo.md)
(which may be brought up by typing
[`?prepareSampleInfo`](https://00berst33.github.io/flowFun/reference/prepareSampleInfo.md)
into the console when `flowFun` has been loaded with
[`library(flowFun)`](https://00berst33.github.io/flowFun/)). Here we
will give a few more thorough examples.

The comparisons must be a named list of named lists. This description
isn’t too helpful on its own, so we will first give a basic example,
then elaborate. The chunk below defines two basic comparisons using the
previously shown table (Example 2).

``` r

comparisons1 <- list(
  old_vs_young = list(age.group = list("<65", "65+")), # comparison 1
  male_vs_female_tmnt = list(treatment = "Treatment", sex = list("female", "male")) # comparison 2
)
```

When defining a named list in R, the left-hand side (LHS) of each `=` is
the name of that item. So, in this example, `comparisons1` has two
elements named `"old_vs_young"` and `"male_vs_female_tmnt"`, each one
corresponding to a test we’d like to make.

Note that these names may be anything you like, but special characters
like `#` and `?` should be avoided. It is also recommended that
underscores `_` or periods `.` are used instead of spaces, to avoid
unexpected behavior from certain functions.

Then for each comparison, we have another named list. This is where we
instruct which groups should be included in that test. In this case,
these lists are `list(age.group = list("<65", "65+"))`, and
`list(treatment = "Treatment", sex = list("female", "male"))`. Here we
must be more careful about list names. The names must also be column
names in your table of sample info, which in this case are `age.group`,
`treatment`, and `sex`. This tells R which columns to check to find the
info it needs to make that particular comparison.

Finally, for each column/variable relevant to that comparison, we
specify which group(s) are of interest. For comparison 1, by setting
`age.group = list("<65", "65+")` we indicate that we’d like to test for
differences between old and young groups (referencing the table above,
this would be sample 1, 5, and 6 vs. 2, 3, and 4). :::{.border
style=“padding: 10px; border: 1px solid \#dee2e6 !important;”}

| sample.name | filename     | age.group |
|:------------|:-------------|:----------|
| Sample 1    | sample_1.fcs | \<65      |
| Sample 2    | sample_2.fcs | 65+       |
| Sample 3    | sample_3.fcs | 65+       |
| Sample 4    | sample_4.fcs | 65+       |
| Sample 5    | sample_5.fcs | \<65      |
| Sample 6    | sample_6.fcs | \<65      |

:::

For comparison 2, by setting `treatment = "Treatment"`, we indicate that
we only want to test within the treatment group, and by setting
`sex = list("female", "male")` we indicate that within that group we
want to test for differences between males and females (this would be
sample 1 and 2 vs. 4 and 6).

| sample.name | filename     | treatment | sex    |
|:------------|:-------------|:----------|:-------|
| Sample 1    | sample_1.fcs | Treatment | female |
| Sample 2    | sample_2.fcs | Treatment | female |
| Sample 4    | sample_4.fcs | Treatment | male   |
| Sample 6    | sample_6.fcs | Treatment | male   |
