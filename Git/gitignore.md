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
