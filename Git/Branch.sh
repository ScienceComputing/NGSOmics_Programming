# Why use branches? 1) Systematically track multiple versions of files; 2) Minimize the risk of conflicting versions; 
# 3) Avoid endless subdirectories; 4) Allow multiple users to work simultaneously

# Display branches
git branch
# *: current branch

# Create a new branch from the current branch
git checkout -b new_branch_name
git branch # Display branches again

# Compare the difference between branches
git diff branch_1 branch_2

# Why switch branches?
# Allow us to work on different components of a project concurrently
git checkout switch_to_another_branch_name
git branch # Confirm if we switch successfully

# Why merge branches?
# The main branch serves as the definitive source of truth for our project. 
# Generally, we allocate each branch to a distinct task. 
# Upon completion of these tasks, it's essential to merge their changes back into the main branch. 
# Doing so ensures that the main branch remains current and accurate.
git merge source_branch_name destination_branch_name
git merge scRNA_seq_preprocess main # Merge the scRNA_seq_preproces branch into the main branch

# Merge conflict
# A conflict arises when a file has differing contents across various branches, making it impossible for them to merge seamlessly into a unified version.
# ***Avoid conflict: use each branch for a specific task; avoid editing a file in multiple branches

# TD: nano
