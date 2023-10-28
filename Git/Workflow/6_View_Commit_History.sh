# Display all commits made to the repo in chronological order, starting with the oldest 
# Show the commit hash, author, date, and commit message
# git log command is very useful for tracking hash and changes at a high level before diving deeper into more specific version control tasks
# Useful to find the commit we want to revert to
git log

# Restrict the number of commits displayed, especially when the project scales
git log -6 # Show the 6 most recent commits

git log -6 file_name # Show the 6 most recent commits of a particular file

# View the commits since 3nd March, 2023
git log --since='Mar 3 2023'

# View the commits between 3nd March, 2023 and 6th March, 2023
git log --since='Mar 3 2023' --until='Mar 6 2023'
