#!/usr/bin/perl -w
use strict;

die "perl $0 [chr] [GeneIDname] > [output.model]" unless @ARGV==2;
my $gff="/vol1/agis/ruanjue_group/wuzhichao/00.data/03.genome/01.oryza/R498/annotationV3/R498_IGDBv3_coreset.gff.msort";
my ($chr,$id)=@ARGV;
my $currentID;

my ($start,$end,$ori);
my %splice;
my $intron_start=0;
open GFF,"$gff" or die $!;
while(<GFF>){
	next if /^#/;
	my @e=split;
	next if $e[0] ne $chr;
	if($e[2] eq "gene"){
		my $tmp_ID=(split /[=;]/,$e[8])[1];
		$tmp_ID=(split /\./,$tmp_ID)[0];
		if ($tmp_ID eq $id){
			$currentID=$id;
			($start,$end)=@e[3,4];
			if($e[6] eq "+"){
				$ori="forward";
			}else{
				$ori="reverse";
			}
		}
	}
	next unless $currentID;
	last if $e[3] > $end+2000;
	if($e[2] eq "CDS" or $e[2]  eq "five_prime_utr" or $e[2] eq "three_prime_utr" or $e[2] eq "exon"){
		my $spliceID=(split /Parent=/,$e[8])[-1];
		my $tmpID=(split /\./,$spliceID)[0];
		next if $tmpID ne $id;
		if ($e[2]  eq "CDS"){
			$e[2]="coding_region";
		}elsif($e[2] eq "five_prime_utr"){
			$e[2]="5' utr";
		}elsif($e[2] eq "three_prime_utr"){
			$e[2]="3' utr";
		}
		$e[3]=$e[3]-$start+1;
		$e[4]=$e[4]-$start+1;
		push @{$splice{$spliceID}{$e[2]}},[@e[3,4]];
		#print "$e[2]\t$e[3] $e[4]\n" if $spliceID eq "OsR498G1018574200.01.T01";
	}
}

my $nsplice=scalar keys %splice;
print "par(mfrow=c($nsplice,1))\n";
my $XAX="T";
foreach my $spliceID (sort keys %splice){
print "#$spliceID\n";
	my $i=0;
	my %exon_start;
	foreach my $pos ( @{$splice{$spliceID}{'exon'}} ){
			my ($a,$b)=@$pos;
			$exon_start{$i}=$a;
			$i++;
	}
	my @exon;my @intron;
	my @i = (sort { $exon_start{$a}<=>$exon_start{$b} } keys %exon_start );
	my $N=scalar @i;
	foreach my $j (0.. $N-1 ){
		my $i=$i[$j];
		my $pos= $splice{$spliceID}{'exon'}[$i];
		my ($a,$b)=@$pos;
		if($j>0){
			my $lastI=$i[$j-1];
			my $intron_start=$splice{$spliceID}{'exon'}[$lastI]->[1] +1;
			my $intron_end=$a-1;
			#print "intron $intron_start $intron_end\n";
			push @intron,"$intron_start-$intron_end";
		}
		#print "exon $a $b\n";
		push @exon,"$a-$b";
	}
	delete $splice{$spliceID}{'exon'};
	my @type_pos;

	print "$spliceID = data.frame( type=c(";
	$i=0;
	foreach my $type ( "5' utr","coding_region","3' utr" ){
		next if(! defined $splice{$spliceID}{$type});
		foreach my $rg ( @{$splice{$spliceID}{$type}} ){
			if($i){
				print  ",\"$type\"";
			}else{
				print  "\"$type\"";
			}
			my $pos =join("-",@$rg );
			push @type_pos,$pos;
			$i++;
		}
	}
	print ",\"exon\"" x scalar @exon;
	print ",\"intron\"" x scalar @intron;
	print "),\n\tcoordinates=c(";
	my @coord=(@type_pos,@exon,@intron);
	(@type_pos,@exon,@intron)=();
	print "\"$_\"," for(@coord[0..$#coord-1]);
	print "\"$coord[-1]\"";
	print ") )\n";
	
	print qq{genemodel.plot(model=$spliceID,start=$start, bpstop=$end, orientation="$ori",xaxis=$XAX)\n\n};
	$XAX="F";
}

