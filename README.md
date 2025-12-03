[![Open in MATLAB Online](https://www.mathworks.com/images/responsive/global/open-in-matlab-online.svg)](https://matlab.mathworks.com/open/github/v1?repo=dnafinder/cuzick&file=cuzick.m)

# cuzick

## 📘 Overview
cuzick is a MATLAB function that implements **Cuzick’s test for trend**, a Wilcoxon-type nonparametric test designed to detect a monotonic trend across three or more **independent groups**.

Unlike multiple pairwise Wilcoxon rank-sum (Mann–Whitney) tests, Cuzick’s test aggregates information across all groups in a single statistic, providing a more powerful and coherent approach when the groups have a natural ordering (e.g., increasing dose levels, ordered risk categories).

The test:

- assumes at least ordinal data,
- requires a meaningful order of the groups,
- uses ranks of the pooled data,
- combines group-wise rank sums with numeric scores that define the trend direction and spacing.

The null hypothesis H₀ is **no trend** across groups, while the alternative hypothesis H₁ is that there is a **monotonic trend** (often increasing) in the distribution/location of the groups, according to the chosen scores.

## ✨ Features
- Simple data structure: `x` is an N-by-2 matrix `[observation, groupLabel]`.
- Group labels must be consecutive integers `1, 2, ..., k` (no gaps), ensuring clear group identification.
- Accepts a vector of **numeric scores** (real values, not necessarily integers) to define:
  - the **ordering** of groups,
  - and the **spacing** between them (e.g., dose levels like 0, 5, 10, 20).
- Computes:
  - Group-wise:
    - number of samples per group,
    - sum of ranks in each group (from pooled ranking),
    - associated scores.
  - Global test quantities:
    - `L` = Σ score(i)·nᵢ
    - `T` = Σ score(i)·Rᵢ (Rᵢ = sum of ranks in group i)
    - `E(T)` and `Var(T)` under H₀
    - z-statistic = |T − E(T)| / √Var(T)
    - one-sided (right-tail) p-value from the standard normal approximation.
- Returns a **structured output** for programmatic use.
- Optional **Display** flag to control printing:
  - `Display = true` (default): prints group summary and test results.
  - `Display = false`: no console output; only returns the STATS struct.

## 📥 Installation
1. Download or clone the repository:
   https://github.com/dnafinder/cuzick

2. Add the folder containing `cuzick.m` to your MATLAB path:
      addpath('path_to_cuzick_folder')

3. Verify that MATLAB can see the function:
      which cuzick

If MATLAB returns the path to `cuzick.m`, the installation is successful.

## ⚙️ Requirements
- MATLAB (any recent version)
- Statistics and Machine Learning Toolbox (required for:
  - `crosstab`
  - `tiedrank`
  - `normcdf`)

## 📈 Usage

Basic call with default integer scores (1, 2, ..., k):

    % x is an N-by-2 matrix:
    %   x(:,1) = observations (ordinal or continuous)
    %   x(:,2) = group labels (1,2,...,k)
    STATS = cuzick(x);

Use custom scores to reflect spacing or coding of the ordered groups:

    % Suppose there are k groups (labels 1..k), but you want scores
    % to reflect actual dose levels, e.g. [0, 5, 10, 20, 40]:
    score = [0 5 10 20 40];
    STATS = cuzick(x, score);

Run the test without printing anything to the Command Window:

    STATS = cuzick(x, [], 'Display', false);

In this mode, the function is “silent” and just returns the STATS structure.

## 🔢 Inputs

cuzick(x)  
cuzick(x, score)  
cuzick(x, score, 'Display', DISPLAY)

- **x**
  - Type: numeric matrix (N×2)
  - Description:
    - `x(:,1)` = observations (at least ordinal, can be continuous)
    - `x(:,2)` = group labels (integer codes)
  - Requirements:
    - Group labels must be **consecutive integers** from 1 to k (1,2,…,k) without gaps.

- **score** (optional)
  - Type: numeric vector (real values), length k
  - Description:
    - Numeric scores associated with each group.
    - May be real-valued (not necessarily integers).
    - Encodes the ordering and spacing of the groups for the trend test.
    - Must not be constant (at least two distinct values).
  - Default:
    - If omitted or empty (`[]`), `score = 1:k`.

- **'Display'** (Name–Value, optional)
  - Type: logical scalar (`true`/`false`) or numeric scalar (`1`/`0`)
  - Default: `true`
  - Description:
    - `true` → prints group table and Cuzick statistics (original behaviour).
    - `false` → no printed output, only the STATS structure is returned.

## 📤 Outputs

cuzick returns a single structured output:

    STATS = cuzick(...);

Fields of **STATS**:

- `STATS.GroupTable`  
  - Type: table  
  - Columns:
    - `Group`     : group label (1, 2, …, k)
    - `Score`     : score assigned to each group (user-defined or default)
    - `Samples`   : number of observations per group (nᵢ)
    - `Ranks_sum` : sum of ranks for each group (Rᵢ), from the pooled ranking of all data

- `STATS.L`  
  - Type: scalar  
  - Description:  
    - L = Σ score(i)·nᵢ (weighted sum of sample sizes)

- `STATS.T`  
  - Type: scalar  
  - Description:  
    - T = Σ score(i)·Rᵢ (weighted sum of rank sums)

- `STATS.E`  
  - Type: scalar  
  - Description:  
    - E(T) under H₀ (no trend)

- `STATS.Var`  
  - Type: scalar  
  - Description:  
    - Var(T) under H₀, including tie correction via `tiedrank`

- `STATS.z`  
  - Type: scalar  
  - Description:  
    - z-statistic = |T − E(T)| / √Var(T), approximated by standard normal under H₀.

- `STATS.pvalue`  
  - Type: scalar  
  - Description:  
    - One-tailed (right-tail) p-value: `p = 1 - normcdf(z)`.

- `STATS.tail`  
  - Type: char/string  
  - Value: `'right'`  
  - Description:
    - Indicates that the test is one-sided and uses the upper tail.

- `STATS.TiesFactor`  
  - Type: scalar  
  - Description:
    - Ties adjustment factor, `2*t`, where `t` is the tie factor returned by `tiedrank`.

## 🧠 Interpretation

- A **large positive z** (and small p-value) suggests evidence for a **monotonic trend** across groups, consistent with the ordering implied by `score`.
- The test is **one-sided (right tail)**:
  - H₀: no trend across groups.
  - H₁: trend in the direction of increasing scores (e.g., increasing severity or dose).
- Small p-values (e.g. p < 0.05) indicate statistically significant evidence of trend.
- The `GroupTable` summarized in STATS helps to inspect:
  - whether higher-scored groups tend to have larger ranks (i.e., larger values),
  - how ties are distributed across groups.

## 📌 Example

Using the mouse metastases example from the help section:

    % Data (metastasis counts)
    d = [0 0 1 1 2 2 4 9 ...
         0 0 5 7 8 11 13 23 25 97 ...
         2 3 6 9 10 11 11 12 21 ...
         0 3 5 6 10 19 56 100 132 ...
         2 4 6 6 6 7 18 39 60];

    % Group labels for cell lines 64, 167, 170, 175, 181 coded as 1..5
    g = [ones(1,8) ...
         2.*ones(1,10) ...
         3.*ones(1,9) ...
         4.*ones(1,9) ...
         5.*ones(1,9)];

    x = [d' g'];

Run Cuzick’s test with de
