# Run the following commands in iTerm2: 
# cd desktop
# perl PerlBasics.pl


$var1=10;
$var2=18;
print "The var1 is: $var1, \n";
print "The var2 is: $var2, \n";
$product=$var1*$var2;
print "The product of var1 and var2 is: $product, \n";
$sum=$var1+$var2;
print "The sum of var1 and var2 is: $sum, \n";
$remainder=$var1%$var2;
print "The remainder of dividing var1 by var2 is: $remainder, \n";

$var1+=9;
print "The new var1 is: $var1, \n";
$var1=10;
$var1-=9;
print "The new var1 is: $var1, \n";
$var1=10;
$var1*=9;
print "The new var1 is: $var1, \n";
$var1=10;
$var1/=9;
print "The new var1 is: $var1, \n";
$var1=10;
$var1=$var1**2;
print "The new var1 is: $var1, \n";

$string="ATC"x3;
print "The sequence is: $string, \n";

# 花括号与条件执行
$abc_x="";
$abc_y=-7;
$abc_x=$abc_y;
if ($abc_x<0) {$abc_x *= -1;}
print "The final value of abc_x is: $abc_x, \n";

# 附加条件
$species = "Homo sapiens";
if ($species eq "Homo sapiens") {
    print "The sequence of origin is human \n";
} elsif ($species eq "Mus musculus") {
    print "The sequence of origin is mouse \n";
} else {
    print "The sequence of origin is neither human nor mouse \n";
}
