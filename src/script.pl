use 5.010;
use POSIX qw(strftime);

sub read_lines {
    my ($file) = @_;

    open my $in, "<:encoding(utf8)", $file or die "$file: $!";
    local $/ = undef;
    my $content = <$in>;
    close $in;
    return split /\n/, $content;
}
my @instances = read_lines("../input/Instances.txt");
my $count = 1;
my $flag = 0;
my $data = strftime "%F", localtime;
open(FIN,">>novasInst.sh");
foreach my $record (@instances) {
	if($record ne "\n"){
		if($flag==0){ #cria a pasta se não houver
                        print FIN ("mkdir ../outputArt/exp-$data-$count\n");
                	$flag=1;
		}
		print FIN ("./projectFinal  -f /home/danilo/codigos/projeto2020/input/instancesGenerated/$record -t 10800 -p 1000 -CC 1 -h 2 -l 0 -m 1 -sp 200000 -ti 0.05677954 -fr 0.99014076 -pv1 0.71070808 -pv2 0.31902837 -pv3 0.02008214 -pv4 0.60150876 -pv5 0.96864417 -tsi 0 >>../outputArt/exp-$data-$count/Resultado-$record.txt 2>&1 \n");	
	}
}

close(FIN);
exit;



