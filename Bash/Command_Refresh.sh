ls -FlhatrR # R shows the nested directories
grep # return contents that match the regex pattern
grep 'lt' file.txt
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
