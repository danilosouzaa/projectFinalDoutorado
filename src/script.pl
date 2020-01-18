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
my @p1 =("0");
my @p2 =("0");
my @p3 = ("1");
my @p4 =("10000");
my @p5 = ("10", "50","100","250", "500", "1000","2500","5000", "7500", "10000");
my @p6 = ("0.1","0.2","0.3","0.4","0.5", "0.6", "0.7", "0.8", "0.9", "-1", "-2");
my $count = 1;
my $flag = 0;
my $data = strftime "%F", localtime;
foreach my $typeSeparation (@p1){
	foreach my $typeLifted (@p2){
		if($typeLifted == 0){
			@p3 =("1");
		}else{
			@p3 = ("1");
		}
		foreach my $minimal (@p3){
			if($typeSeparation==1){
                        	foreach my $record (@instances) {
                                	if($record ne "\n"){
										if($flag==0){
                                       	        	system("mkdir ../output/exp-$data-$count");
                                       	        	$flag=1;
										}
                                  	        open(FIN,">>../output/exp-$data-$count/Resultado-$record.txt");
        									print FIN ("Guloso exp. $record minimal $minimal tipoLifted $typeLifted  \n");
                                        	close(FIN);
                                        	system("./projectFinal instancesGenerated/$record 0 1000 -CC $typeSeparation $typeLifted $minimal 0 0 0 >>../output/exp-$data-$count/Resultado-$record.txt 2>&1");
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
                                        	                                system("mkdir ../output/exp-$data-$count");
                                        	                                $flag=1;
                                        	                        }
                                        	                        open(FIN,">>../output/exp-$data-$count/Resultado-$record.txt");
                                        	                        print FIN ("Grasp exp. $record minimal $minimal tipoLifted $typeLifted szPool $pool nIteration $iteGrasp alpha $alpha \n");
                                        	                        close(FIN);
                                        	                        system("./projectFinal instancesGenerated/$record 0 1000 -CC $typeSeparation $typeLifted $minimal $pool $iteGrasp $alpha >>../output/exp-$data-$count/Resultado-$record.txt 2>&1");
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

