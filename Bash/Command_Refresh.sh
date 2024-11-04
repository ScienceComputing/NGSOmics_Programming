ls -FlhatrR 
# -F adds helpful markers to each item in the list: a “/” after a name indicates it’s a directory, an asterisk “*” shows it’s a program, and if there’s no symbol, it’s a regular text file.
# -l prints long format information about files including file mode, link count, owner, group, size, date, time, name
# -h print
# -a lists all files and directories, including the hidden ones that begin with a dot
# t sorts the output by modification time, with the most recently modified files appearing first
# r reverses the sorting order, so the oldest files appear first
# -R shows the nested directories

grep # return contents that match the regex pattern
grep 'lt' file.txt
cd 1_dir/2_dir/ && cat file.csv | grep 'abc'
grep '[ltde]' file.txt # [] constitues a matching set - e.g., return anything containing l/t/d/e; ^[] constitutes an inverse set


cat # concatenate contents line by line
tail -n # return only the last -n lines
head -n # return only the first -n lines
wc -w # return the word count
wc -l # return the line count
sed # replace string based on the matched pattern
regex # regular expressions used to filter files, data within files, match arguments
# test the regex using regex101.com
sort | uniq -c # sort the contents alphabetically, and then do the unique count; uniq only considers adjacent for its distinct operation 
cat file.txt | sort | uniq -c | head -n 6 # return the top 6 most frequent objects 
