# Compare an unstaged file with the last committed version
git diff file_name

# Compare a staged (git add) file with the last committed version
git diff -r HEAD file_name

# Compare multiple staged files with the last committed versions
git diff -r HEAD

# Display the differences between the current state of the specified file and its state in the previous commit (the commit just before the current one)
# Helpful for reviewing changes made to the file in the most recent commit
git diff -r HEAD~1 file_name

# Display the differences between the current state of the specified file and its state two commits ago
# Useful for seeing how a particular file has changed over the course of two commits
git diff -r HEAD~2 file_name
