## Why use gitignore
.gitignore is particularly useful if we are running many programs, where the number of files tracked is fairly large.

## Commonly ignored files
- APIs
- Credentials
- System files
- Software dependencies

## Go to the root of our local Git folder

```
touch .gitignore
vim .gitignore
```

## Insite the .gitignore file, write the files or directories we want Git not to track 
```
scRNA_draft.*
```
The pattern `scRNA_draft.*` matches any file or directory whose name begins with `scRNA_draft.`.


```
/scRNA_draft.*
```
The pattern `scRNA_draft.*` matches any file or directory whose name begins with `scRNA_draft.`, such as `scRNA_draft.qmd`, but not `code/scRNA_draft.qmd`.


```
/scRNA
```
The pattern `scRNA` will match a directory `scRNA` and paths underneath it, but will not match a regular file or a symbolic link  for `scRNA`.


```
/scRNA/code_experiment
scRNA/code_experiment
```
The leading slash is not relevant if there is already a middle slash in the pattern.


```
scRNA/*
```
The pattern `scRNA/*`, matches `scRNA/matrix.mtx.gz` (a regular file), `scRNA/code` (a directory), but it does not match `scRNA/code/estimate_rna_velocity.py` (a regular file), as the asterisk in the pattern does not match `code/estimate_rna_velocity.py` which has a slash in it.
