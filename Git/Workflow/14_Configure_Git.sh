# Access a list of customizable settings
git config --list



# Settings for one specific project
git config --local

# Settings for all projects
git config --global
git config --global user.email XXX@gmail.com
git config --global user.name 'YYY UUU'

# Settings for every user on this computer
git config --system

# !Create an alias for committing files by executing ci
config --global alias.ci 'commit -m'
# Commit files 
git ci F100.fastq.gz 'Upload the FASTQ file'

# !Create an alias for unstaging files by executing unstage
git config --global alias.unstage 'reset HEAD'
# The alias should not overwrite the existing git or shell command

# Show all the created aliases
git config --global --list
