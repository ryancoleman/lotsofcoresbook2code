#!/usr/bin/perl
#

$[ = 0;                 # set array base to 1
$, = ' ';               # set output field separator
$\ = "\n";              # set output record separator

#die "Usage: fold_prop.pl\n" unless @#ARGV == 0;

while (<STDIN>)
{
  chop;
  @Fld = split(' ', $_);
    
  $tag = $Fld[2] if ($Fld[0] eq "TAG");
  $tag =~ s/\"//g;
  $fn  = $Fld[2] if ($Fld[0] eq "FUNCTION");
  $fn  =~ s/\"operator/\"/;

  if ($Fld[0] eq "EXPR")
  {
#    printf "%s %s\n", $tag, $fn;
    print STDOUT << "EOF1";
// $tag
template <>
struct TagVisitor<$tag, PrintTag> : public ParenPrinter<$tag>
{ 
  static void visit($tag op, PrintTag t) 
    { t.os_m << $fn; }
};
EOF1
  }
}

  
  
