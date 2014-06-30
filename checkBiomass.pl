use strict;
use XML::LibXML;
use Data::Dumper;

sub parseGroupCpds{

	my $g = shift;
	my $h = {};
	for my $spr ($g->getElementsByTagName("speciesReference")){
		$h->{$spr->getAttribute("species")} = $spr->getAttribute("stoichiometry");
	}
	return $h;

}

sub parseModel{

	my $f = shift;
	my $parser = XML::LibXML->new;
	my $doc = $parser->parse_file($f);
	my $cH = {};
	my $rH = {};

	# compounds
	for my $cpd ($doc->getElementsByTagName("species")){
		my $id = $cpd->getAttribute("id");
		$cH->{$id}{id}=$id;
		$cH->{$id}{name}=$cpd->getAttribute("name");
		$cH->{$id}{compartment}=$cpd->getAttribute("compartment");
		for my $p ($cpd->getElementsByTagName("p")){
			$p->textContent =~ /(\w+)\: +(\S+)/;
			$cH->{$id}{lc($1)}=$2;
		}
		$cH->{$id}{elements} = parseElement($cH->{$id}{formula});
	}

	# reactions
	for my $rxn ($doc->getElementsByTagName("reaction")){
		my $id = $rxn->getAttribute("id");
		$rH->{$id}{id} = $id;
		$rH->{$id}{name} = $rxn->getAttribute("name");
		$rH->{$id}{reversible} = $rxn->getAttribute("reversible");

		for my $group ($rxn->getElementsByTagName("listOfReactants")){
			$rH->{$id}{reactants} = parseGroupCpds($group);
		}
		for my $group ($rxn->getElementsByTagName("listOfProducts")){
			$rH->{$id}{products} = parseGroupCpds($group);
		}
	}

	return ($cH,$rH);

}


sub parseElement{

	my $f = shift;
	my @elementList = ('C','H','O','N','P','S','K','Na','Cl','Ca','Mg','Mn','Co','Cu','Se','Fe');
	my $elementCount = {};
	
	for my $e (@elementList){
		if ($f =~ m/($e)[A-Z]/){
			$elementCount->{$e} = 1;
		}elsif ($f =~ m/($e)$/){
			$elementCount->{$e} = 1;
		}elsif ($f =~ m/$e(\d+)n/){
			$elementCount->{$e} = $1*3;
		}elsif ($f =~ m/$e(\d+)/){
			$elementCount->{$e} = $1;
		}else{
			$elementCount->{$e} = 0;
		}
	}
	
	return $elementCount;
}


sub checkMass{

	my ($r, $cH) = @_;
	my $reactantsElements = {};
	my $productsElements = {};
	my $elementsAll = {};
	my $f = 0;
	for my $reac (sort keys %{$r->{reactants}}){
		for my $e (keys %{$cH->{$reac}{elements}}){
			$reactantsElements->{$e} += $cH->{$reac}{elements}{$e} * $r->{reactants}{$reac};
			$elementsAll->{$e}++;
		}
	}
	for my $prod (sort keys %{$r->{products}}){
		for my $e (keys %{$cH->{$prod}{elements}}){
			$productsElements->{$e} += $cH->{$prod}{elements}{$e} * $r->{products}{$prod};
			$elementsAll->{$e}++;
		}
	}

	for my $e (sort keys %$elementsAll){
		$elementsAll->{$e} = $reactantsElements->{$e} - $productsElements->{$e};
		$f = 1 if $elementsAll->{$e} != 0;
	}

	return ($f, $elementsAll);

}

sub checkCharge{

	my ($r, $cH) = @_;
	my $reactantsCharge = 0;
	my $productsCharge = 0;
	my $chargeAll = 0;
	my $k = 0;
	for my $reac (sort keys %{$r->{reactants}}){
		$reactantsCharge += $cH->{$reac}{charge};
	}
	for my $prod (sort keys %{$r->{products}}){
		$productsCharge += $cH->{$prod}{charge};
	}

	
	$chargeAll = $reactantsCharge - $productsCharge;
	$k = 1 if $chargeAll != 0;
	

	return ($k, $chargeAll);

}

sub checMassChargeConsistency{
	my ($cH, $rH) = @_;
	my $incMassRH = {};
	my $incChargeRH = {};
	for my $rxn (sort keys %$rH){
		if ($rxn !~ /^R_EX_/){
			my ($f, $diffAtoms) = checkMass($rH->{$rxn}, $cH);
			my ($k, $diffCharge) = checkCharge($rH->{$rxn}, $cH);

			$incMassRH->{$rxn}{diff} = $diffAtoms if ($f);
			$incChargeRH->{$rxn}{diff} = $diffCharge if ($k);
		}

	}
	return ($incMassRH, $incChargeRH);
}

sub printRxnFormula{
	my ($r, $cH) = @_;
	print $r->{id},": ",$r->{name},"\treversible: ",$r->{reversible},"\n";
	my @leftSide; my @rightSide;
	my @leftSideFormulae; my @rightSideFormulae;
	for my $c (sort keys %{$r->{reactants}}){
		push @leftSide, ($r->{reactants}{$c} == 1) ? $c : $r->{reactants}{$c}.' '.$c;
		push @leftSideFormulae, ($r->{reactants}{$c} == 1) ? $cH->{$c}{formula} : $r->{reactants}{$c}.' '.$cH->{$c}{formula};
	}
	for my $c (sort keys %{$r->{products}}){
		push @rightSide, ($r->{products}{$c} == 1) ? $c : $r->{products}{$c}.' '.$c;
		push @rightSideFormulae, ($r->{products}{$c} == 1) ? $cH->{$c}{formula} : $r->{products}{$c}.' '.$cH->{$c}{formula};
	}
	if ($r->{reversible} eq 'false'){
		print join(" + ", @leftSide), " => ", join(" + ", @rightSide),"\n";
		print join(" + ", @leftSideFormulae), " => ", join(" + ", @rightSideFormulae),"\n";
	}elsif ($r->{reversible} eq 'true'){
		print join(" + ", @leftSide), " <=> ", join(" + ", @rightSide),"\n";
		print join(" + ", @leftSideFormulae), " <=> ", join(" + ", @rightSideFormulae),"\n";
	}

	return 1;
}

sub printIncMassRxns{
	my ($incR, $r, $c) = @_;
	my $count = 0;
	for my $rid (sort keys %$incR){
		$count++;
		print "inconsistent Mass $count: ";
		printRxnFormula($r->{$rid}, $c);
		for my $e (sort keys %{$incR->{$rid}{diff}}){
			print $e,": ",$incR->{$rid}{diff}{$e},"\n" if $incR->{$rid}{diff}{$e};
		}
	breakLine();

	}
}

sub breakLine{
	my $line = '-' x 50;
	print $line,"\n";

}

my $f = 'Ca_iYZ766_30-Jun-2014.xml';
my ($cHash, $rHash) = parseModel($f);
my ($incMassRxnsHash, $incChargeRxnsHash) = checMassChargeConsistency($cHash, $rHash);
printIncMassRxns($incMassRxnsHash, $rHash, $cHash);
# print Dumper(parseElement('CNO'));

