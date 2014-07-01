use strict;


sub parseElement{

	my $f = shift;
	my @elementList = ('C','H','O','N','P','S','K','Na','Cl','Ca','Mg','Mn','Co','Cu','Se','Fe','Zn');
	my $elementCount = {};
	
	for my $e (@elementList){
		if ($f =~ m/($e)[A-Z]/){
			$elementCount->{$e} = 1;
		}elsif ($f =~ m/($e)$/){
			$elementCount->{$e} = 1;
		}elsif ($f =~ m/$e(\d+)n/){
			$elementCount->{$e} = $1*3;
		}elsif ($f =~ m/$e(\d+\.\d+)/){
			$elementCount->{$e} = $1;
		}elsif ($f =~ m/$e(\d+)/){
			$elementCount->{$e} = $1;
		}else{
			$elementCount->{$e} = 0;
		}
	}
	
	return $elementCount;
}



sub calMW{
	my $elementNum = shift;
	my $elementMW = {
	'C'=> 12.0107,
	'H'=> 1.00794,
	'O'=> 15.9994,
	'N'=> 14.0067,
	'P'=> 30.9738,
	'S'=> 32.066
	};
	my $mw = 0;
	for my $e (sort keys %$elementMW){
		$mw += $elementMW->{$e}*$elementNum->{$e};
	}
	return $mw;
}


sub printFormula{
	my $elementNum = shift;
	my @atoms = qw/C H O N P S/;
	
	for (@atoms){
		print $_,int(10000*$elementNum->{$_})/10000 if (int(10000*$elementNum->{$_})/10000);
	}
	return 1;
}

sub breakLine{
	my $line = '-' x 50;
	print $line,"\n";

}

my $biomassElementNum = {C => 36.96841758,
H => 31.7000695500001,
N => 7.67926869999999,
O => 20.1416073700001,
P => 1.01065100000001,
S => 0.9763001,
};
print "Biomass molecular weight = ",calMW($biomassElementNum),"\n";
print "Biomass formula: "; printFormula($biomassElementNum); print "\n";
breakLine();


my $constituentsFormulae = {
'M_carbo_met_c' => 'C37.038H61.73O30.865',
'M_dna_met_c' => 'C31.879H36.74O19.442N11.834P3.242',
'M_peptido_met_c' => 'C40H62O21N8',
'M_plipid_met_c' => 'C55.0129H102.8475O13.029N0.397P1.415',
'M_protein_met_c' => 'C43.855H12.806O21.273N10.746S1.999',
'M_rna_met_c' => 'C29.7H33.21O21.589N11.778P3.107',
'M_teich_met_c' => 'C158.5581H373.333O251.9457N3P48.1860',
'M_trace_met_c' => 'C29.06H34.914O16.931N9.379P2.28S0.199',
};


my $constituentsPercent = {
'M_carbo_met_c' => 4.32,
'M_dna_met_c' => 2.6,
'M_peptido_met_c' => 10.09,
'M_plipid_met_c' => 7.6,
'M_protein_met_c' => 52.84,
'M_rna_met_c' => 6.55,
'M_teich_met_c' => 8,
'M_trace_met_c' => 4.94,
};


for my $c (sort keys %$constituentsFormulae){
	my $eCount = parseElement($constituentsFormulae->{$c});
	print "$c molecular weight = ",my $mw = calMW($eCount),"\n";
	print "$c formula: "; printFormula($eCount); print "\n";
	print "$c original % = ", $constituentsPercent->{$c}/100,"\tnow = ",
	 int(10000*$constituentsPercent->{$c}/100/$mw*1000)/10000,"\n";
	breakLine();
}