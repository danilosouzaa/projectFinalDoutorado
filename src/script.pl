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
my @p1 =("0","1");# 0- Grasp, 1 - Greedy
my @p2 =("0","1"); # 0 - letchford, 1 -Ballas
my @p3 = ("0", "1"); # 1 - minimal, 0- no minimal
my @p4 =("200000");## poolSize
my @p5 =("200000"); #number of iteration
my @p6 = ("-1");#alpha
my $count = 1;
my $flag = 0;
my $data = strftime "%F", localtime;
foreach my $typeSeparation (@p1){
	foreach my $typeLifted (@p2){
		if($typeLifted == 0){
			@p3 =("1");
		}else{
			@p3 = ("0","1");
		}
		foreach my $minimal (@p3){
			if($typeSeparation==1){
                        	foreach my $record (@instances) {
                                	if($record ne "\n"){
						if($flag==0){ #cria a pasta se não houver
                                       	        	system("mkdir ../outputArt/exp-$data-$count");
                                       	        	$flag=1;
						}
                                  	        open(FIN,">>../outputArt/exp-$data-$count/Resultado-$record.txt");
        					print FIN ("Guloso exp. $record minimal $minimal tipoLifted $typeLifted  \n");
                                        	close(FIN);
                                        	system("./projectFinal benchmark/$record 10800 1000 -CC $typeSeparation $typeLifted $minimal 200000 200000 -1 >>../outputArt/exp-$data-$count/Resultado-$record.txt 2>&1");
					}
				}
				$count++;
                        	$flag = 0;
			}else{
				foreach my $pool (@p4){
					foreach my $iteGrasp (@p5){
						foreach my $alpha (@p6){
                               				foreach my $record (@instances) {
                               					if($record ne "\n"){
                                        				if($flag==0){
                                        	                                system("mkdir ../outputArt/exp-$data-$count");
                                        	                                $flag=1;
                                        	                        }
                                        	                        open(FIN,">>../outputArt/exp-$data-$count/Resultado-$record.txt");
                                        	                        print FIN ("Grasp exp. $record minimal $minimal tipoLifted $typeLifted szPool $pool nIteration $iteGrasp alpha $alpha \n");
                                        	                        close(FIN);
                                        	                        system("./projectFinal benchmark/$record 10800 1000 -CC $typeSeparation $typeLifted $minimal $pool $iteGrasp $alpha >>../outputArt/exp-$data-$count/Resultado-$record.txt 2>&1");
                                        	                }
                                        	     	}
                                        	        $count++;
                                        	        $flag = 0;
                                        	}
					}
                                }
                        }
		}
	}
}
exit;

