#script that implements naive version of DBSCAN to explore clustering of stems
#some functions borrowed from detect_regions_stem
#currently not forming stems til after clustering
#input: Rfam family name
#output: consensus stems
#run from Cluster-graph/consensus-scripts

use strict;
use Common;
use POSIX;

if (scalar(@ARGV) < 1) {
    die "Please input a Rfam family name, along with the number of seqs you want (optional)\n";
} else {
    my $sfold_profile = 0;
    my $print_heatmaps = 0;
    my $totalseq = 0; #total number of seqs we want, use 0 for all in process_alignmt
    my $totalseq = $ARGV[1] if scalar(@ARGV)>1; 
    my $aln_file = "../../data/analysis/$ARGV[0]"."_analysis/process_alignmt.out";
    my $nat = process_native($aln_file,$totalseq);
    my $type = "Sfold";
    if ($sfold_profile) {$type = run_sfold_profiling($ARGV[0],$nat)};
#    my $lengths = process_alignmt($aln_file);
#    my $stems = process_stems($ARGV[0],$lengths);
    my $newclus = 0;
    my $feats;
    my $repeat = 1;
    my $found = {};
    while ($repeat) {
	my $data = process_profiling($ARGV[0],$nat,$type);
	$feats = $data->[0];
	my $lengths = $data->[1];
	my $cells = process_trip_feat($feats,$lengths,$found);
	scatterplot_cells($cells,$ARGV[0]."-init") if $print_heatmaps;
	my $eps = find_eps($lengths);
	my $clusters = cluster($cells,$eps);
#	scatterplot_cluster($clusters->[0],$ARGV[0]);
	my $clusstems = print_clusters($clusters->[0],$cells,$feats,$nat,$ARGV[0]);
#    my $clusstems = print_clusters($clusters,$cells,$stems,$nat);
	$newclus = find_centroid($clusstems->[0],$lengths,$ARGV[0],$data->[2]);
#	print_clus_stats($newclus);
	improve_cluster($newclus,$data,$ARGV[0],$nat,$clusstems->[1],$clusters->[1]);
#	print_clus_stats($newclus);
	if (found_newclus($found,$newclus,$feats)) {
	    $found = get_truepos($newclus,$feats);
#	    print_clus_stats($newclus);
	    $type = resample($found,$ARGV[0],$feats);
	    $type = "Resample";
	    print "Resampling:\n";
	} else {
	    $repeat = 0;
	}
    }
    print "FINAL CLUSTERS:\n";
    print_clus_stats($newclus);
    print_clus($newclus);
    print_clus_to_file($newclus,$ARGV[0],$totalseq);
#    my $rfam = get_Rfam_ref($ARGV[0]);
#    native_stems($rfam);
}

sub print_clus_stats {
    my $newclus = shift;
    for (my $i = 1; $i<scalar(@$newclus); $i++) {
	next if (!$newclus->[$i]);
	my $clus = $newclus->[$i];
	my $feats = $clus->{"constraints"};
	printf "found %d seqs for cluster %d\n", scalar(keys %$feats),$i;
    }
}

#grabs the native helices
sub process_native {
    my ($file,$totalseq) = @_;
    my $seq;
    my $fam;
    my %nat;
    my $thresh = 0.2; #accept if below this thresh
    open ALN, "<$file" or die "cannot open $file\n";
    foreach (<ALN>) {
	if (/Sequence \.\.\/seq\/(\S+)\/(\S+)\.txt/) {
	    $seq = $2;
	    $fam = $1;
	} elsif (/Profile ([\d\s\!]+)/) {
#	    print "processing $seq $1\n";
	    my @nat = split(/\s/,$1);
	    $nat{$seq} = \@nat;
	}
    }
    close(ALN);
    $totalseq = 0 if ($totalseq < 0 || $totalseq > scalar(keys %nat));
    if ($totalseq) {
	my @keys = keys %nat;
	my @passed = ();
	while (scalar(keys %nat) > $totalseq) {
	    my $rand = rand();
#	    print "Found random number $rand\n";
	    my $del = shift(@keys);
	    if ($rand > $thresh) {
		print "Deleting $del\n";
		delete $nat{$del};
	    } else {
		push(@passed,$del);
	    }
	    if (scalar(@keys)==0) {
		@keys = @passed;
		@passed = ();
	    }
	}
    }
    print "Total seqs: ",scalar(keys %nat),"\n";
    return \%nat;
}

#runs profiling on sfold outputs of fam
sub run_sfold_profiling {
    my ($fam,$nat) = @_;
    my $run_sfold = 1;
    foreach my $seq (keys %$nat) {
	my $seqfile = "../../seq/$fam/$seq".".txt";
	my $profile = "../../data/profiling_output/Rfam_$fam/$seq/$seq"."_sfold.out";
	my $sfold = "../../data/sfold_output_Rfam/$fam"."/$seq"."/sample_1000.out";
#run sfold
	my $outdir = "/scratch/tmp/$seq" . "_output";

	if ($run_sfold) {
	    print "Running Sfold on $seq\n";
	    `/usr/bin/sfold -o $outdir $seqfile`;
	    my $outfile = $outdir."/sample_1000.out";
	    $profile = "../../data/profiling_output/Rfam_$fam/$seq";
	    `mkdir $profile`;
	    $profile = $profile . "/$seq"."_sfold.out";
	    $sfold = "../../data/sfold_output_Rfam/$fam"."/$seq";
	    `mkdir $sfold`;
	    $sfold = $sfold."/sample_1000.out";
	    `mv $outfile $sfold`;
	}
	
	`../../RNAStructProfiling/RNAprofile -sfold $sfold -g -v $seqfile > $profile`;
    }
    return "Sfold";
}

#runs consensus_DP and grabs lengths
sub process_alignmt {
    my $file = shift;
    my %lengths;
    my @norm;
    my $seq;
    my $ave;
    my $med;
    my @out = `perl consensus_DP.pl $file`;
    foreach (@out) {
	if (/Length\: (\d+)/) {
	    $lengths{$seq} = $1;
#	    print "processing $seq\n";
	} elsif (/(\S+)\:/) {
	    $seq = $1;
#	} elsif (/max length is ([\d\.]+), ([\d\.]+), (\d+), (\d+)/) {
#	    $med = $1;
#	    $ave = $2;
#	    $lengths{"min"} = $3;
#	    $lengths{"max"} = $4;
#	    $lengths{"range"} = $4-$3;
	}
    }
#    foreach (values %lengths) {
#	my $val = abs($med-$_)*($med/$_);
#	push(@norm,$val);
#    }
#    $lengths{"median"} = $med;
#    $lengths{"stdev"} = stdev(\@norm);
#    $lengths{"stdev"} = stdev([values %lengths]);
#    print "Stdev of lengths is ",$lengths{"stdev"},"\n";
    return \%lengths;
}

sub process_stems {
    my ($fam,$nat) = @_;
    my @seqs = keys %{$nat};
    my %all;
    foreach my $seq (@seqs) {
	my $file = "../../data/profiling_output/Rfam_$fam"."/$seq".".out";
	my @out = `perl calc_stem_class.pl $file`;
	my %stems;
#	print "\nprocessing $seq ";
	foreach (@out) {
	    if (/Stem (\d+) has ([\d+\s]+) \(([\d+\s]+)/) {
		$stems{$1} = [$2,$3];
#		print "$1 ";
	    }
	}
	$all{$seq} = \%stems;
    }
    return \%all;
}

#grabs feature triplet info and puts into hash
#processes input profiling file based on type:
# 0 = gtboltzmann input, 1 = sfold input, 2 = resample input
sub process_profiling {
    my ($fam,$nat,$type) = @_;
    my @seqs = keys %{$nat};
    my %feats;
    my %lengths;
    my %nucl;
    my %files;
    foreach my $seq (@seqs) {
	my @feat = ();
	my $file = "../../data/profiling_output/Rfam_$fam"."/$seq".".out";
	if ($type eq "Sfold") {
	    $file = "../../data/profiling_output/Rfam_$fam"."/$seq"."/$seq"."_sfold.out";
	} elsif ($type eq "Resample") {
	    $file = "../../data/profiling_output/Rfam_$fam/$seq/resample.out";
	    my $found = `test -e $file && echo Found || echo Missing`;
	    if ($found =~ /Missing/) {
		$file = "../../data/profiling_output/Rfam_$fam"."/$seq"."/$seq"."_sfold.out";
	    }
	}
	open OUT, "<$file" or die "cannot open $file\n";
	foreach (<OUT>) {
	    if (/Featured helix (\d+)\: (\d+ \d+ \d+) with freq (\d+)/) {
#	    $feats{$1} = $2;
		$feat[$1] = $2;
	    } elsif (/seq in .+ is ([UCAG]+) with length (\d+)/) {
		$nucl{$seq} = [0,split(//,$1)];
		$lengths{$seq} = $2;
	    }
	}
	close(OUT);
	$feats{$seq} = \@feat;
	$files{$seq} = $file;
    }
    return [\%feats,\%lengths,\%nucl,\%files];
}

#takes all stems from all seqs, and normalizes HC trip info
#takes normalized triplets, integerizes them and makes pairs hash
#cell hash returned to be used with dbscan to group
#return featcells to use features, cells to use stems
sub process_trip {
    my ($stems,$feats,$lengths) = @_;
    my @seqs = keys %$lengths;
    my $median = median([values %$lengths]);
    my %cells;
    my %featcells;
    my @heatmap = (0) x $median;
    foreach my $seq (@seqs) {
	my $mystems = $stems->{$seq};
	my $length = $lengths->{$seq};
	my $myfeats = $feats->{$seq};
	foreach my $stem (keys %$mystems) {
	    my @hcs = split(/\s/,$mystems->{$stem}[0]);
	    my $label = "$seq $stem";
	    foreach my $feat (@hcs) {
		my $featlabel = "$seq $feat";
		my @trip = split(/\s/,$myfeats->[$feat]);	
		my $i = ceil($trip[0]*$median/$length);
		my $j = ceil($trip[1]*$median/$length);
#		print "processing $label $feat @trip -> ($i,$j)\n";
		for (my $k = 0; $k < $trip[2]; $k++) {
		    add(\%cells,"$i $j",$label);
		    add(\%featcells,"$i $j",$featlabel);
		    $i++;
		    $j--;
		}
	    }
#	    $stocells{$label} = \@mycells;
	}
#	last;
    }
#    print "Found seqs for cells ",keys %cells,"\n";
#    return \%cells;
    return \%featcells;
}

#same as process_trip but instead of using stems, use features
#compare against already found features in %found; add if not present
sub process_trip_feat {
my ($feats,$lengths,$found) = @_;
    my $median = median([values %$lengths]);
    my %featcells;
#    my @heatmap = (0) x $median;
    foreach my $seq (keys %$feats) {
	my $myfeats = $feats->{$seq};
	my $length = $lengths->{$seq};
	my @coords = ();
	if (exists $found->{$seq}) {
	    @coords = @{$found->{$seq}};
	}
	for (my $feat = 1; $feat<scalar(@$myfeats);$feat++) {
	    my $featlabel = "$seq $feat";
	    my $coord = $myfeats->[$feat];
	    remove_char(\@coords,$coord);
	    add_feature($coord,$median,$length,\%featcells,$featlabel);
	}
#add those features from before that might not have made it through resampling
	if (scalar(@coords)) {
	    my $k = scalar(@$myfeats);
	    foreach my $coord (@coords) { 
		add_feature($coord,$median,$length,\%featcells, "$seq $k");
		push(@$myfeats,$coord);
		print "Adding $seq ($k) $coord\n";
		$k++;
	    }
	}
    }
#    print "Found seqs for cells ",keys %cells,"\n";
    return \%featcells;
}

#takes coordinates and adds to featcells
sub add_feature {
    my ($coord,$median,$length,$featcells,$featlabel) = @_;
    my @trip = split(/\s/,$coord);	
    my $i = ceil($trip[0]*$median/$length);
    my $j = ceil($trip[1]*$median/$length);
#		print "processing $label $feat @trip -> ($i,$j)\n";
    for (my $k = 0; $k < $trip[2]; $k++) {
	add($featcells,"$i $j",$featlabel);
	$i++;
	$j--;
    }
}

#finds an appropriate radius eps to use in dbscan
#based on density of lengths; need to adjust for normalization
#add median to lengths
sub find_eps {
    my $lengths = shift;
    my $median = median([values %$lengths]);
    my $n = scalar(keys %$lengths);
    my $eps = 0;
    my $minpts = scalar(keys %$lengths)/4;
    my %vals = map {$_ => ceil(abs($median-$lengths->{$_})*$median/$lengths->{$_})} keys %$lengths;
    my %ltoseq;
    foreach (keys %$lengths) {
	my $coord = ceil(abs($median-$lengths->{$_})*$median/$lengths->{$_});
	add(\%ltoseq,$coord,$_);
    }
    my $noise = $n;
    my $label; #label{length} = cluster
    while ($noise > $n/2) {
	$label = DBscan(\%ltoseq,$eps,$minpts,"length");
	my @noise = map {$label->{$_} == 0 ? ($_):() } keys %$label;
	my @n = map {scalar(@{$ltoseq{$_}}) } @noise;
	$noise = sum(\@n);
	print "found eps $eps with $noise seq as noise\n";
	$eps++;
    }
    print_eps($label,\%ltoseq);
    $lengths->{"median"} = $median;
    $eps = 2*($eps-1); #preparing to be used for 2D, so double
    return [$eps,$minpts];
}

sub print_eps {
    my ($label,$ltoseq) = @_;
    my %clus;
    foreach my $l (keys %$label) {
	add(\%clus,$label->{$l},$l);
    }
    foreach my $k (keys %clus) {
	my $lengths = $clus{$k};
#	print "For cluster $k: ";
	foreach my $len (@$lengths) {
	    my $seqs = $ltoseq->{$len};
#	    print "@$seqs ($len), ";
	}
    }
}

#clusters the cells, given the eps found from length densities
#if it produces over half of the points as noise, lower minpts
sub cluster {
    my ($cells,$params) = @_;
    my $eps = $params->[0];
    my $minpts = $params->[1];
    my $n = scalar(keys %$cells);
    my $noise = $n;
    my $label;
    while ($noise > $n/2) {
	if ($minpts < 3) {
	    print "minpts below critical level of 3; raising eps to ",$eps+1,"\n";
	    $minpts = $params->[1];
	    $eps = 3;
	}
	$label = DBscan($cells,$eps,$minpts,"coord");
	my @noise = map {$label->{$_} == 0 ? ($_):() } keys %$label;
	my @n = map {scalar(@{$cells->{$_}}) } @noise;
	$noise = sum(\@n);
	print "found minpts $minpts with $noise points as noise out of $n\n";
	$minpts--;
    }
    $minpts++;
    return [$label,$minpts];
}


#$ltoseq{pt} = [members], $label{pt} = cluster
sub DBscan {
    my ($ltoseq,$eps,$minpts,$type,$maxradius) = @_;
    my $C = 0;
    my %label;
    my @sorted = sort {scalar(@{$ltoseq->{$b}}) <=> scalar(@{$ltoseq->{$a}}) } keys %$ltoseq;
    foreach my $p (@sorted) {
	#print "considering seq $_\n";
	next if (exists $label{$p});
	my $neighbors;
	$neighbors = rangeQuery($p,$eps,$ltoseq) if ($type eq "length");
	$neighbors = rangeQuery2($p,$eps,$ltoseq) if ($type eq "coord");
	if ($neighbors->[1] < $minpts) {
	    $label{$p} = 0;
	    next;
	}
	$C++;
	$label{$p} = $C;
	my @set = @{$neighbors->[0]};
	next if (!scalar(@set)); #if all neighbors are same point, no extra neighbors to process
	foreach my $pt (@set) {
#	    print "checking point $pt\n";
	    if (exists $label{$pt}) {
		$label{$pt} = $C if (!$label{$pt});
		next;
	    }
	    $label{$pt} = $C;
	    $neighbors = rangeQuery($pt,$eps,$ltoseq) if ($type eq "length");
	    $neighbors = rangeQuery2($pt,$eps,$ltoseq) if ($type eq "coord");
	    push(@set,@{$neighbors->[0]}) if ($neighbors->[1] >= $minpts);
	}
    }
#    print "Total num clusters: $C\n";
    return \%label;
}

#finds all lengths within eps of val, adds to neighbors
#one dimensional
sub rangeQuery {
    my ($val,$eps,$ltoseq) = @_;
    my $numneigh = scalar(@{$ltoseq->{$val}});
    my @neighbors = ();
    for (my $i = 1; $i <= $eps; $i++) {
	my $v = $val+$i;
	if (exists $ltoseq->{$v}) {
	    $numneigh += scalar(@{$ltoseq->{$v}});
	    push(@neighbors,$v);
	}
	$v = $val - $i;
	if (exists $ltoseq->{$v}) {
	    $numneigh += scalar(@{$ltoseq->{$v}});
	    push(@neighbors,$v);
	}
    }
#    print "found for $val $numneigh neighbors: @neighbors\n";
    return [\@neighbors,$numneigh];
}

#version of rangeQuery that's two dimensional
#also, we count the number of seqs represented not just points
sub rangeQuery2 {
    my ($val,$eps,$ltoseq) = @_;
    my %found = map {(split(/\s/,$_))[0] => 1} @{$ltoseq->{$val}};
#    my $numneigh = scalar(@{$ltoseq->{$val}});
    my @neighbors = ();
    my @val = split(/\s/,$val);
    for (my $x = $eps; $x >= $eps*-1; $x--) {
	my $i = $val[0]+$x;
	for (my $y = $eps-$x; $y >= $x-$eps; $y--) {
	    next if (!$x && !$y); #already processed original cell
	    my $j = $val[1]+$y;
	    my $v = "$i $j";
	    if (exists $ltoseq->{$v}) {
#		$numneigh += scalar(@{$ltoseq->{$v}});
		map {$found{(split(/\s/,$_))[0]}++} @{$ltoseq->{$v}};
		push(@neighbors,$v);
	    }
	}
    }
    my $numneigh = scalar(keys %found);
#    print "found for $val $numneigh neighbors: @neighbors\n";
    return [\@neighbors,$numneigh];
}

#takes the labels and prints them for matlab scatter plot
sub scatterplot_cluster {
    my ($label,$fam) = @_;
    my $file = "../../Matlab/data/dbscan_plot_$fam".".txt";
    open OUT,">$file" or die "cannot open $file\n";
    foreach (keys %$label) {
	print OUT "$_ $label->{$_}\n";
    }
    close(OUT);
    return;
}

sub scatterplot_cells {
    my ($cells,$fam) = @_;
    my $file = "../../Matlab/data/dbscan_cells_$fam".".txt";
    open OUT,">$file" or die "cannot open $file\n";
    foreach (keys %$cells) {
	print OUT "$_ ",scalar(@{$cells->{$_}}),"\n";
    }
    close(OUT);
    return;
}

#prints the stems making up each cluster
#returns an array indexed by cluster num, containing the stems in each cluster
sub print_clusters {
    my ($label,$cells,$stems,$nat,$fam) = @_;
    my %clus;
    my %found;
    my $verbose = 0;
    my $print_heatmap = 1;
    if ($print_heatmap) {
	my $file = "../../Matlab/data/dbscan_init_$fam".".txt";
	open OUT,">$file" or die "cannot open $file\n";
    }
    my $feat = 1; #if we want to cluster feats, instead of stems
    foreach my $cell (keys %$label) { #gives info by every cluster
	add(\%clus,$label->{$cell},$cell);
    }
    my @clusters;
    my @cluscells;
    foreach my $clust (sort {$a <=> $b} keys %clus) {
	print "For cluster $clust: \n" if ($verbose);
	my %allstems;
	foreach my $cell (@{$clus{$clust}}) {
	    my $stms = $cells->{$cell};
	    map {$allstems{$_}++} @$stms;
	    if ($print_heatmap && $clust) {
		print OUT "$cell ",scalar(@$stms),"\n";
	    }
	}
	my @sorted = sort {$a cmp $b} keys %allstems;
	my %myclus;
	foreach my $stem (@sorted) { #for every stem in the cluster
	    my @stem = split(/\s/,$stem);
	    my $native = $nat->{$stem[0]};
	    my $steminfo;
	    my $steminfoprint;
	    if ($feat) { #if using features to cluster
		my $coords = $stems->{$stem[0]}[$stem[1]];
		if (exists_num($native,$stem[1])) {
		    $steminfoprint = "(*$stem[1]) [$coords]";
		    $steminfo = "$stem[1] | $coords";
		} else {
		    $steminfoprint = "($stem[1]) [$coords]";
		    $steminfo = "$stem[1] | $coords";
		}
	    } else { #if using stems to cluster
		my $coords = $stems->{$stem[0]}{$stem[1]};
		$steminfo = "$coords->[0] | $coords->[1]";
		my @hcs = map {exists_num($native,$_) ? ("$_*") : ($_)} (split(/\s/,$coords->[0]));
		$steminfoprint =  "(@hcs) [$coords->[1]]";
	    }
	    print "\t$stem[0] $steminfoprint\n" if ($verbose);
	    add(\%myclus,$stem[0],$steminfo);
	    $found{$stem} = $clust if ($clust);
	}
#	map {$found{$_} = $clust} @sorted if ($clust);
	$clusters[$clust] = \%myclus;
    }
    close OUT if ($print_heatmap);
    return [\@clusters,\%found];
}


#combines all stems of a seq into a superstem
#finds the median coordinates of all seq stems
#variation: use normalized coords?
sub find_centroid {
    my ($clusters,$lengths,$fam,$nucls) = @_;
    my @newclus;
    my %hctoclus;
    my $median = $lengths->{"median"};
    my $verbose = 1;
    for (my $clus = 0; $clus<scalar(@$clusters); $clus++) {
	#skip if its the the noise cluster = 0
	if (!$clus) {
	    $newclus[$clus] = $clusters->[$clus];
	    next;
	}
	print "For cluster $clus:\n" if ($verbose);
	my $myclus = $clusters->[$clus];
	my @sorted = sort {$a cmp $b} keys %$myclus;
	my %newclus;
	my %seqs;
	my %constrain;
	my @I;
	my @J;
	my @K;
	my @L;
	my @lengths;
	my @longest;
	my @mfe;
	my @bpnums;
	foreach my $seq (@sorted) {
	    my $stems = $myclus->{$seq};
	    my $nucl = $nucls->{$seq};
	    my @is = ();
	    my @js = ();
	    my $longest = [0,0,0];
	    my @feats = ();
	    my $bpnum = 0;
	    my $minHC = 0;
	    my $mincoords;
	    foreach my $stem (@$stems) {
		my @stem = split(/ \| /,$stem);
		my @coords = split(/\s/,$stem[1]);
		if ($stem[0] < $minHC || !$minHC) {
		    $minHC = $stem[0];
		    $mincoords = \@coords;
		}
		my $k = scalar(@coords)-1; #allowing for both feat (i,j,k) and stem (i,j,k,l) coords
		push(@is,($coords[0],$coords[0]+$coords[2]-1)); #push the boundaries, to calc k and l
		push(@js,($coords[1],$coords[1]-$coords[$k]+1));
		$bpnum += $coords[2]; #assuming using feats (i,j,k)
		if ($longest->[2] < $coords[2]) {$longest = [@coords];}
		elsif ($longest->[2] == $coords[2]) { push(@$longest,@coords);}
		push(@feats,$stem[0]);
		
	    }
	    my ($i,$j,$k,$l) = @{find_stem_coords(\@is,\@js)};
	    push(@I,$i);
	    push(@J,$j);
	    push(@K,$k);
	    push(@L,$l);
	    push(@longest,$longest);
	    push(@lengths,$lengths->{$seq});
	    my $name = "$fam"."_c$clus" . "_$seq";
	    my $bounds =  [$i,$i+$k-1,$j-$l+1,$j];
#	    my $mfe = run_mfe($name,$seq,$nucl,$bounds);
#	    push(@mfe,$mfe->[1]);
	    push(@bpnums,$bpnum);
	    $seqs{$seq} = [$i,$j,$k,$l];
	    $constrain{$seq} = \@feats;
	    print "\t$seq (@feats) [$i,$j,$k,$l] {$bpnum bp}\n" if ($verbose); #can print (L: @$longest) or $mfe->[1]
	}
	$newclus{"seqs"} = \%seqs;
	$newclus{"constraints"} = \%constrain;
	my @ks = map {$_->[2]} @longest;
	my $centroid = [median(\@I),median(\@J),median(\@K),median(\@L)];
	$newclus{"centroid"} = $centroid;
	$newclus{"offset"} = median(\@lengths);
#calculating coordinates for putative median
	my $offset = $median-$newclus{"offset"};
	my $medcoords = find_medcoords($centroid,$median,$newclus{"offset"});
	$newclus{"medcoords"} = $medcoords;
	$newclus{"bpnum"} = median(\@bpnums);
	print "   Median offset: $offset\n" if ($verbose);
	print "   Median longest HC: ",median(\@ks),"\n" if ($verbose);
	print "   Median centroid: @$centroid\n" if ($verbose); 
	print "   Median bpnum: ",$newclus{"bpnum"},"\n" if ($verbose);
	print "   Med coords: @$medcoords\n" if ($verbose);
#	print "   Median MFE: ",median(\@mfe),"\n" if ($verbose);
#print out the missing seqs
	my @seqs = keys %$lengths;
	remove_char(\@seqs,"median");
	my @foundseqs = keys %seqs;
	my @missing = map {exists_char(\@foundseqs,$_) ? () : ($_)} @seqs;
	print "   Missing seqs:\n\t@missing\n" if ($verbose);
	$newclus{"missing"} = \@missing;
	$newclus[$clus] = \%newclus;
#	last if ($clus == 3);
    }
    return \@newclus;
}

sub find_stem_coords {
    my ($is,$js) = @_;
    my $i = min($is);
    my $j = maxm($js);
    my $k = maxm($is) - $i + 1;
    my $l = $j - min($js) + 1;
    return [$i,$j,$k,$l];
}

#finds what the putative centroid of median length would be
#shift to putative centroid by rate offset/clus_med
sub find_medcoords {
    my ($centroid,$median,$clus_med) = @_;
    my $offset = $median-$clus_med;
    my $upperi = ceil($centroid->[0]+ $centroid->[0]*($offset/$clus_med));
    if ($upperi < 1) { $upperi = 1;}
    my $loweri = $upperi +$centroid->[2]-1;
    my $lowerj = floor($centroid->[1]+$centroid->[1]*($offset/$clus_med));
    if ($lowerj > $median) {$lowerj = $median;    }
    my $upperj = $lowerj - $centroid->[2]+1;
    return [$upperi,$loweri,$upperj,$lowerj];
#    my $med_coords = [$centroid->[0]+$offset,$centroid->[0]+$offset+$centroid->[2]-1,$centroid->[1]+$offset-$centroid->[3]+1,$centroid->[1]+$offset];
}

#based on centroid, go look for missing seq hc that might fit the bill
#look within range of 0 to offset from theoretical median
#offset composed of normalized length difference, plus length of hc to allow for split hc
sub improve_cluster {
    my ($newclus,$data,$fam,$native,$found,$minpts) = @_;
    my $feats = $data->[0];
    my $lengths = $data->[1];
    my $median = $lengths->{"median"};
    my @seqs = keys %$lengths;
    remove_char(\@seqs,"median");
    my %combine;
    for (my $i = 1; $i<scalar(@$newclus); $i++) {
	my $clus = $newclus->[$i];
	my $size = scalar(keys %{$clus->{"seqs"}});
	print "For Cluster $i\n";
	if ($size < $minpts) {
	    my $badclus = $newclus->[$i];
	    $newclus->[$i] = 0;
	    print "\tCluster $i with ",$size," below $minpts threshold: deleting\n";
	    next;
	}
	my ($score,$newconstraints) = @{find_missing(@_,$i)};
	my $total = $score + $size; #since each present member of clus gets score of 1, so sum = size
	if ($total <= 0) {
	    $newclus->[$i] = 0;
	    print "\tCluster $i with score of $total: deleting\n";
	} else {
	    #add found seq HC to cluster: update constraints, 
	    my $constraints = $clus->{"constraints"};
	    my $stemcoords = $clus->{"seqs"};
	    foreach my $seq (keys %$newconstraints) {
		my $found = $newconstraints->{$seq};
		$constraints->{$seq} = $found;
		my $myfeats = $feats->{$seq};
		my @is = ();
		my @js = ();
		foreach my $feat (@$found) {
		    my $max = $myfeats->[$feat];
		    my @coords = split(/\s+/,$max);
		    push(@is,($coords[0],$coords[0]+$coords[2]-1)); #push the boundaries, to calc k and l
		    push(@js,($coords[1],$coords[1]-$coords[2]+1));
		}
		$stemcoords->{$seq} = find_stem_coords(\@is,\@js);
		print "\tAdding $seq (@{$constraints->{$seq}}) [@{$stemcoords->{$seq}}]\n";
	    }
	    my $newcentroid = recalc_centroid($stemcoords);
	    $clus->{"centroid"} = $newcentroid;
	    print "\tNew centroid: @$newcentroid\n";
	    print "\tTotal score of cluster $i: $total\n";
	}
#	last if ($i == 3);
    }
}

#recalculates the centroid after missing seqs added in
sub recalc_centroid {
    my $stemcoords = shift;
    my @is;
    my @js;
    my @ks;
    my @ls;
    foreach my $seq (keys %$stemcoords) {
	my $coords = $stemcoords->{$seq};
	push(@is,$coords->[0]);
	push(@js,$coords->[1]);
	push(@ks,$coords->[2]);
	push(@ls,$coords->[3]);
    }
    return [median(\@is),median(\@js),median(\@ks),median(\@ls)];
}

sub find_mfe {
    my ($med_coords,$lengths,$median,$nucls,$name,$seqs) = @_;
    my %mfe;
    foreach my $seq (@$seqs) {
	my $length = $lengths->{$seq};
	my $diff = $length-$median;
	my $nucl = $nucls->{$seq};
	my $bounds = find_boundaries($med_coords,$diff,$length);
	my $mfe = run_mfe($name."_$seq",$seq,$nucl,$bounds);
	$mfe{$seq} = $mfe;
	my $bp = 0;
	map {$_ eq '('? $bp++ : ()} @{$mfe->[0]};
	print "\tFor $seq:\t",@{$mfe->[0]}," [$mfe->[1]] ($bp bp)\n";
    }
    return \%mfe;
}

# found{seq} = clus
# clus{seq} = [i,j,k,l]
#return in %constrain list of features found
sub find_missing {
    my ($newclus,$data,$fam,$native,$found,$minpts,$i) = @_;
    my $feats = $data->[0];
    my $lengths = $data->[1];
    my $clus = $newclus->[$i];
    my $median = $lengths->{"median"};
    my $missing = $clus->{"missing"};
    my $centroid = $clus->{"centroid"};
    my $bpnums = $clus->{"bpnum"};
#    my $clusmed = $clus->{"offset"};
#    my $offset = $median-$clusmed;
#    my $med_coords = 
#	[$centroid->[0]+$offset,$centroid->[0]+$offset+$centroid->[2]-1,$centroid->[1]+$offset-$centroid->[3]+1,$centroid->[1]+$offset];
    my $med_coords = $clus->{"medcoords"};
    my $verbose = 1;
    my $totalscore = 0;
    my %constrain;
    print "For cluster $i:\n" if ($verbose);
    foreach my $seq (@$missing) {
	print "   For seq $seq, " if ($verbose);
	my $feat = $feats->{$seq};
	my $length = $lengths->{$seq};
	my $nat = $native->{$seq};
	my $diff = $length-$median;
	my $nucl = $data->[2]{$seq};
	my $bounds = find_boundaries($med_coords,$diff,$length);
	my $scores = search_window($seq,$bounds,$nucl,$centroid,$found,$nat,$diff,$fam,$newclus,$i,$lengths,$data->[3]);
#	shrink_window($bounds,$scores,$median);
#	choose_best(\@scores,$found,$seq,$i,$newclus,$verbose,$nat);
	my $mfe = run_mfe($fam,$i,$seq,$nucl,$bounds);
	my $bpnum = 0;
	if ($mfe->[1] != 0) {
	    foreach my $bp (@{$mfe->[0]}) {
		if ($bp eq '(') {$bpnum++;}
	    }
	    if ($bpnum >= $bpnums) {$totalscore++;}
#	    elsif ($bpnum >= ceil($bpnums/2)) {$scores{$seq} = 0;}
	    elsif ($bpnum < ceil($bpnums/2)) {$totalscore--;}
	} else {
#	    $scores{$seq} = -2;
	    $totalscore -= 2;
	}

#	my $adjusted = translate_mfe($mfe,$bounds);
#	my $mincoords = match_to_hc($adjusted,$scores);
#	$mincoords = get_constraints($adjusted) if (!$mincoords);
#	$constrain{$seq} = $mincoords if ($mincoords);
#	print @{$mfe->[0]}," [$mfe->[1]] {$bpnum bp}\n" if ($verbose);
#	print "\t$seq: ";

#if the found is a feature, save it
	my @found = map {$_->[0] < scalar(@$feat) ? ($_->[0]) : ()} @$scores;
	if (scalar(@found)<1) {next;}
	$constrain{$seq} = \@found;
	print "found potential features @found\n" if ($verbose);
    }
    print "Cumulative score of missing seqs: $totalscore\n";
    return [$totalscore,\%constrain];
}

#right now finding only the innermost HC to constrain;
#can potentially constrain entire MFE
#returns ref to array of (i,j,k)
sub get_constraints {
    my $adjusted = shift;
    my $lastbp = 0;
    my $k = 1;
    foreach my $bp (@$adjusted) {
	if ($lastbp) {
	    if ($bp->[0] == $lastbp->[0]-1 && $bp->[1] == $lastbp->[1]+1) {
		$k++;
		$lastbp = $bp;
	    } else {
		push(@$lastbp,$k);
#		print "found a break at @$bp so using @$lastbp\n";
		return $lastbp;
	    }
	} else {
	    $lastbp = $bp;
	}
    }
    return 0 if (!$lastbp);
    push(@$lastbp,$k);
    return $lastbp;
}

#bounds = [i i' j' j]
sub run_mfe {
    my ($fam,$i,$seq,$nucl,$bounds) = @_;
    my $verbose = 0;
    my $mkdirs = 0;
    my $name = "$fam"."_c$i" . "_$seq";    
#make sequence file to run MFE
    my $dir = "../../consensus/seqs/";
    $dir = mkdirs($dir,$fam,$i) if ($mkdirs);
    my $file = $dir . $name . ".txt";
    open SEQ, ">$file" or die "cannot open $file\n";
    print SEQ ">$name, $seq for window @$bounds\n";
    my @copy = @$nucl;
    my $length;
    my $diff = $bounds->[2]-$bounds->[1];
    my $k = $bounds->[1]-$bounds->[0]+1;
    if ($diff < 4) {
	$length = $bounds->[3]-$bounds->[0]+1;
	my @part = splice(@copy,$bounds->[0],$length);
	print SEQ @part;
    } else {
	my @part = splice(@copy,$bounds->[0],$k);
	print SEQ @part;
	print SEQ "\nAAAA\n";
	@copy = @$nucl;
	my $l = $bounds->[3]-$bounds->[2]+1;
	@part = splice(@copy,$bounds->[2],$l);
	print SEQ @part;
	$length = $k+$l+4;
    }
    close(SEQ);
#make constraints file
    $dir = "../../consensus/constraints/";
    $dir = mkdirs($dir,$fam,$i) if $mkdirs == 1;
    my $constraints = $dir . $name . ".txt";
    make_constraints($constraints,$diff,$k,$length);
    my $gtmfe = "../../Desktop/gtfold-master/bin/gtmfe";
    my $paramdir = "../../Desktop/gtfold-master/gtfold-mfe/data/Turner99";
    my @mfe = `$gtmfe --paramdir $paramdir -c $constraints $file`;
    print "@mfe\n" if $verbose;
    my ($bp,$mfe);
    foreach (@mfe) {
	if (/^[\.\(\)]+$/) {
	    chomp;
	    my @bp = split(//);
	    if ($diff > 3) {
		splice(@bp,$k,4,'x','x','x','x');
	    } elsif ($diff > 0) {
		my @x = 'x' x ($diff-1);
		splice(@bp,$k,$diff-1,@x);
	    } 
	    $bp = \@bp;
	} elsif (/Minimum Free Energy\:\s+(\-?[\d\.]+)/) {
	    chomp;
	    $mfe = $1;
	}
    }
    my $ct = "$name" . ".ct";
    `mv $ct ../../consensus/stems/`;
    return [$bp,$mfe];
}

#makes dirs under $dir for $fam, then cluster $i
sub mkdirs {
    my ($dir,$fam,$i) = @_;
    $dir = $dir . "$fam/";
    `mkdir $dir`;
    $dir = $dir . "cluster$i/";
    `mkdir $dir`;
    return $dir;
}

#calculate bp, then translate back to original coords 
sub translate_mfe {
    my ($mfe,$bounds) = @_;
    my $bp = $mfe->[0];
    my @left;
    my @bps;
    for (my $i = 0; $i < scalar(@$bp); $i++) {
	if ($bp->[$i] eq '(') {
	    push(@left,$i+1);
	} elsif ($bp->[$i] eq ')') {
	    my $left = pop(@left);
	    push(@bps,[$left,$i+1]);
	}
    }
    my $diff = $bounds->[2]-$bounds->[1];
    my $i_offset = $bounds->[0]-1;
    my $j_offset = $bounds->[0]-1;
    if ($diff > 3) {
	$j_offset = $bounds->[2]-($bounds->[1]-$bounds->[0]+1)-5;
    }
    my @adjusted = map {[$_->[0]+$i_offset,$_->[1]+$j_offset]} @bps;
    return \@adjusted;
}

#matches the mfe bp to found HC
sub match_to_hc {
    my ($adjusted,$scores) = @_;
    my %mfebp;
    map {$mfebp{"@$_"}++} @$adjusted;
    my $found=0;
    my @mincoords;
    foreach my $hc (@$scores) {
	for (my $k = 0; $k<$hc->[3]; $k++) {
	    my $i = $hc->[1] + $k;
	    my $j = $hc->[2] - $k;
	    if (exists $mfebp{"$i $j"}) {
		@mincoords = @$hc if (!$found); 
		$found = $hc->[0];
		print "\tFound most freq HC $found in mfe at @$hc\n";
		last;
	    }
	}
    }
    return 0 if (!$found);
    return [splice(@mincoords,1,3)];
}

#forbid pairing within wings, and with dummy filler inbetween wings
sub make_constraints {
    my ($constraints,$diff,$k,$length) = @_;
    open CON,">$constraints" or die "cannot open $constraints\n";
    my $overlap = 0;
#prohibit any filler
    if ($diff > 3) {
	print CON "P ",$k+1," 0 4\n";
    } elsif ($diff > 0) {
	print CON "P ",$k+1," 0 ",$diff-1,"\n" if ($diff != 1);
    } else {
	$overlap = abs($diff)+1;
    }
#prohibit diagonals involving first nuc
    for (my $j = $k-$overlap; ($j-3)/2 >= 1; $j--) {
	my $n = floor(($j-3)/2);
	print CON "P 1 $j $n\n";
    }
#prohibit all other diagonals, starting with other nucleotides
    for (my $i = 2;$i<$k-$overlap-3;$i++) {
	my $n = floor(($k-$overlap-$i+1-3)/2);
	print CON "P $i ",$k-$overlap," $n\n";
    }
#determine where 3' wing begins
    my $jp;
    if ($diff > 3) {
	$jp = $k + 5;
    } elsif ($diff > 0) {
	$jp = $k + $diff;
    } else {
	$jp = $k+1;
    }
#do same thing for other wing
    for (my $j = $length; ($j-$jp+1-3)/2 >= 1; $j--) {
	my $n = floor(($j-$jp+1-3)/2);
	print CON "P $jp $j $n\n";
    }
    for (my $i = $jp+1; $i<$length-3; $i++) {
	my $n = floor(($length-$i+1-3)/2);
	print CON "P $i $length $n\n";
    }
    close(CON);
}

#searchs the window within bounds by making dotplot
#make sure you get the profiling file right
sub search_window {
    my ($seq,$bounds,$nucl,$centroid,$found,$nat,$diff,$fam,$allclus,$i,$lengths,$files) = @_;
    my @scores;
    my $verbose = 0;
    my $move = 0; #to move from another cluster; disabled for now
    my $newclus = $allclus->[$i];
    my $med_coords = $newclus->{"medcoords"};
    my $length = $lengths->{$seq};
    my $median = $lengths->{"median"};
    my $file = $files->{$seq};
#    my $dotplot = make_dotplot($bounds,$nucl);
    print "For seq $seq ($diff), looking within bounds [@$bounds]\n" if ($verbose);
#    my $file = "../../data/profiling_output/Rfam_$fam"."/$seq".".out";
    open OUT,"<$file" or die "cannot open $file\n";
    foreach (<OUT>) {
	chomp;
	if (/Helix (\d+) is (\d+) (\d+) (\d+) .+ with freq (\d+)/) {
	    my $i_end = $2+$4-1;
	    my $j_end = $3-$4+1;
	    if ($bounds->[0]<= $2 && $i_end <= $bounds->[1] && $bounds->[2]<= $j_end && $3 <= $bounds->[3]) {
		my $coords = [$2,$3,$4,$4];
		my $score;
		my $ratios = calc_single_poisson($median,$med_coords,$coords,$length);
		my $coverage = ($4/$centroid->[2] + $4/$centroid->[3])/2;
#		    my $avepoisson = ave($ratios);
		printf "\tFound potential $_ with poisson %.2f (%.2f %.2f %.2f) coverage %.2f",
		ave($ratios),$ratios->[0],$ratios->[1],$ratios->[2],$coverage if ($verbose);
#		printf "\tFound potential $_ with coverage %.2f",$coverage if ($verbose);
		if ($move) {
		    if (exists $found->{"$seq $1"}) { #if in another cluster, pick the larger one
			my $old = $found->{"$seq $1"};
			print " ++ ($old)" if ($verbose);
			my $oldclus = $allclus->[$old];
			if (scalar(keys %{$newclus->{"seqs"}}) > scalar(keys %{$oldclus->{"seqs"}})) {
			    $found->{"$seq $1"} = $i;
			    $newclus->{"seqs"}->{$seq} = $oldclus->{"seqs"}->{$seq};
			    delete $oldclus->{"seqs"}->{$seq};
			    print " (moving)" if ($verbose);
			}
		    }
		}
		if ($verbose) {
		    print " ** " if (exists_num($nat,$1));
#		    print "\n";
		}
#		    if ($coverage > 1) {
#			$score = ave($ratios)+$5/1000 + 1;
#		    } else {
#			$score = ave($ratios)+$5/1000 + $coverage;
#		    }
		if (implausible_HC($ratios,$coords,$med_coords,$length,$median)) {
		    print " xxx\n" if ($verbose);
		    next;
		}
		print "\n" if ($verbose);
		my $info = [$1,$2,$3,$4,$5,@$ratios,$coverage];
		push (@scores,$info);
#		    push(@scores,[$score,$1,$coords]);
	    }
	}
    }
#    print_dotplot($dotplot,$bounds);
#    return $dotplot;
    return \@scores;
}

#allowed opposite offset of 5 or 10% of offset
#implausible if exceeds any allotted offset
#or if all of poisson ratios below 0.5
#ignore if the offset is less than twice the length of the 
sub implausible_HC {
    my ($ratios,$coords,$med_coords,$length,$median) = @_;
#    return 1 if (maxm($ratios) < 0.5);
    my $k = $med_coords->[1]-$med_coords->[0];
    my $l = $med_coords->[3]-$med_coords->[2];
    my $offset = $length-$median;
#if a small offset, can trust poisson
    if (abs($offset) < $k || abs($offset) < $l) {
	return 1 if (maxm($ratios) < 0.5);
	return 0;
    }
#if a large offset, need to calculate actual indels
    my $verbose = 0;
    print "\nCoords @$coords, med @$med_coords " if ($verbose);
    my ($insert,$del,$myinsert,$mydel) = (0,0,0,0);
    my $opp_offset = ceil($offset/10);
    if (abs($opp_offset) < 5) {$opp_offset = 5;}
    if ($offset < 0) {
	$del = $offset - $opp_offset;
	$insert = $opp_offset;
    }
    else {
	$insert = $offset + $opp_offset;
	$del = -1*$opp_offset;
    }
#the required number of indels to align i's
    my $pre_i = $coords->[0]-$med_coords->[0];
    print "before is $pre_i, " if ($verbose);
    if ($pre_i < 0) { $mydel = $pre_i;}
    else {$myinsert = $pre_i;}
#required number of indels between i and j
    my $between = (($coords->[1]-$coords->[0])-($med_coords->[3]-$med_coords->[0]));
    print " between is $between, " if ($verbose);
    if ($between < 0) { $mydel += $between;}
    else {$myinsert += $between;}
#required number indels after j
    my $after = (($length-$coords->[1])-($median-$med_coords->[3]));
    print " after is $after" if ($verbose);
    if ($after < 0) {$mydel += $after;}
    else {$myinsert += $after;}
#have we exceeded indel allotment?
    if (($myinsert - $insert) > 0) {
	print "\n Too many insert: $myinsert, $insert" if ($verbose); 
	return 1;}
    if (($mydel - $del) < 0) {
	print "\n Too many del: $mydel, $del" if ($verbose); 
	return 1;}
    return 0;
}

sub match {
    my ($i,$j) = @_;
    my $match = 0;
    if ($i eq 'A' || $i eq 'a') {
	$match = 1 if ($j eq 'U' || $j eq 'u' || $j eq 'T' || $j eq 't');
    } elsif ($i eq 'C' || $i eq 'c') {
	$match = 1 if ($j eq 'G' || $j eq 'g');
    } elsif ($i eq 'U' || $i eq 'u' || $i eq 'T' || $i eq 't') {
	$match = 1 if ($j eq 'G' || $j eq 'g' || $j eq 'A' || $j eq 'a');
    } elsif ($i eq 'G' || $i eq 'g') {
	$match = 1 if ($j eq 'U' || $j eq 'u' || $j eq 'T' || $j eq 't' || $j eq 'C' || $j eq 'c');
    }
    return $match;
}

#picks the best candidate stem
#if belongs to another cluster already, place in larger cluster
#change alliance
#if moving from old cluster, inactivate it in old one so we can erase deleted clusters
#if the missing seqs we've found belong to a bigger, mark so if majority of missings belong elsewhere, delete current cluster
sub choose_best {
    my ($scores,$found,$seq,$i,$newclus,$verbose,$nat) = @_;
    my $clus = $newclus->[$i];
    my @sorted = sort {$b->[0] <=> $a->[0]} @$scores;
    my $best;
    while (!exists $clus->{$seq} && scalar(@sorted)) {
	$best = shift(@sorted);
#	print "Considering @$best\n" if ($verbose);
	if (exists $found->{"$seq $best->[1]"}) {
	    my $old = $found->{"$seq $best->[1]"};
	    my $oldclus = $newclus->[$old];
	    if (scalar(keys %$clus) > scalar(keys %$oldclus)) {
		$clus->{$seq} = $best->[2];
		$found->{"$seq $best->[1]"} = $i;
		delete $oldclus->{$seq};
		$oldclus->{"deleted"}++;
		print "Changing $seq $best->[1] from clus $old to $i\n" if ($verbose);
	    }
	} else {
	    $clus->{$seq} = $best->[2];
	}
    }
    if ($verbose) {
	if (exists $clus->{$seq}) {
	    print "   Chose HC $best->[1] with score $best->[0]";
	    print " **" if (exists_num($nat,$best->[1]));
	    print "\n";
	} else {
	    print "No viable candidate\n";
	    $best = 0;
	}
    }
    return $best;
}

#finds the boundaries to look for missing seq's HCs
#boundaries for the beginning of the stem to the end of the stem
sub find_boundaries {
    my ($medcentroid,$diff,$length) = @_;
    my ($upperi,$loweri,$upperj,$lowerj);
    my $range = 0.25;
    my $pseudo = 2;
    if ($diff < 0) {
	$upperi = floor($medcentroid->[0]+($diff*(1+$range)));
	$loweri = ceil($medcentroid->[1]-($diff*$range));
	$upperj = floor($medcentroid->[2]+($diff*(1+$range)));
	$lowerj = ceil($medcentroid->[3]-($diff*$range));
    } elsif ($diff > 0) {
	$upperi = floor($medcentroid->[0]-($diff*$range));
	$loweri = ceil($medcentroid->[1]+($diff*(1+$range)));
	$upperj = floor($medcentroid->[2]-($diff*$range));
	$lowerj = ceil($medcentroid->[3]+($diff*(1+$range)));
    } else {
	$upperi = $medcentroid->[0]-$pseudo;
	$loweri = $medcentroid->[1]+$pseudo;
	$upperj = $medcentroid->[2]-$pseudo;
	$lowerj = $medcentroid->[3]+$pseudo;
    }
    $upperi = 1 if ($upperi < 1);
    $upperj = 1 if ($upperj < 1);
    $lowerj = $length if ($lowerj > $length);
    return [$upperi,$loweri,$upperj,$lowerj];
}

#aggregates all the true positive features found in the cluster together
#store as the max triplets;
sub get_truepos {
    my ($newclus,$feats) = @_;
    my %found;
    for (my $i = 1; $i<scalar(@$newclus);$i++) {
	next if (!$newclus->[$i]);
	my $clus = $newclus->[$i];
	my $constraints = $clus->{"constraints"};
	foreach my $seq (keys %$constraints) {
	    my $found = $constraints->{$seq};
	    my $feat = $feats->{$seq};
	    #use $found if just want feat num, not coords
	    if (exists $found{$seq}) { 
		my $myfound = $found{$seq};
		my @coords = map {exists_char($myfound,$feat->[$_]) ? ():($feat->[$_])} @$found;
		push(@$myfound,@coords);
	    }
	    else {
		my @coords = map {$feat->[$_]} @$found;
		$found{$seq} = \@coords;}
	}
    }
    return \%found;
}

#if one of the new clusters contains all new features never part of a cluster before...
#don't worry about counting the nonfeatures found in a missing seq search; the cluster's feats will disqualify
# every single cluster has to have at least one previously found feat for this to be a NO
#new find stored in newclus, found contains the previously found 
sub found_newclus {
    my ($found,$newclus,$feats) = @_;
    return 1 if (scalar(keys %$found) == 0);
    my $verbose = 1;
    for (my $i = 1; $i<scalar(@$newclus); $i++) {
	next if (!$newclus->[$i]);
	my $clus = $newclus->[$i];
	my $constrain = $clus->{"constraints"};
	my $prior = 0;
	foreach my $seq (keys %$constrain) {
	    my $coords = $found->{$seq};
	    my $feat = $feats->{$seq};
	    my $clusfeat = $constrain->{$seq};
#	    print "Checking @$clusfeat ($feat->[$clusfeat->[0]]) against @$coords\n";
	    my @exists = map {exists_char($coords,$feat->[$_]) ? ($feat->[$_]) : ()} @$clusfeat;
	    if (scalar(@exists)) {
		print "Found for overlap for $seq, cluster $i not new\n" if $verbose;
		$prior = 1;
		last;
	    }
	}
	if (!$prior) {
	    print "Cluster $i is new: resample\n" if $verbose;
	    return 1;
	}
    }
    print "No new clusters found\n" if $verbose;
    return 0;
}
 
#forbid the features shown to be false positives, ie not part of remaining clusters
sub resample {
    my ($found,$fam,$feats) = @_;
    my %found;
    my $mkdir = 1;
#prohibit all features not found in clusters
    foreach my $seq (keys %$feats) {
	my $feat = $feats->{$seq};
	my $myfound = [];
	$myfound = $found->{$seq} if (exists $found->{$seq});
	print "examining $seq with found @$myfound\n";
	if (scalar(@$feat) <= scalar(@$myfound)+1) {next;}
	my $name = "$fam" . "_$seq";
	my $cfile = "../../consensus/constraints/$fam/$name" . ".txt";
	open CON, ">$cfile" or die "cannot open $cfile\n";	
	for (my $i = 1; $i<scalar(@$feat); $i++) {
	    if (exists_char($myfound,$feat->[$i])) { next;}
	    print CON "P $feat->[$i]\n";
	}
	close(CON);
#run sfold with constraints
	my $seqfile = "../../seq/$fam/$seq".".txt";
	my $outdir = "/scratch/tmp/$name" . "_output";
	`/usr/bin/sfold -o $outdir -f $cfile $seqfile`;
#test new structures are energetically plausible
	my $output = "../../data/gtboltzmann_outputs/Rfam_$fam/$seq/output.samples";
	$output = "../../data/sfold_output_Rfam/$fam"."/$seq"."/sample_1000.out";
	my $stats = energy_calcs($output);
	my $out = $outdir . "/sample_1000.out";
	my $newstats = energy_calcs($out);
	if (!plausible($stats,$newstats)) {
	    print "implausible new sample for $seq; sticking with original\n";
	    next;
	}
#run profiling on new sfold output
	my $profile = "../../data/profiling_output/Rfam_$fam/$seq/resample.out";
	`../../RNAStructProfiling/RNAprofile -sfold $out -g -v $seqfile > $profile`;
#save files
	my $sfolddir = "../../data/sfold_output_Rfam/$fam/$seq/resample/";
	`mkdir $sfolddir` if $mkdir;
	`mv $out $sfolddir`;
    }
    return "Resample";
}

sub energy_calcs {
    my ($out) = @_;
    my @fe;
    my $sfold = 1;
    if ($out =~ /gtboltzmann/) {$sfold = 0;}
    open OUT,"<$out" or die "cannot open $out\n";
    if ($sfold) {    
	foreach (<OUT>) {
	    if (/\s\d+\s+(\-[\d\.]+) /) {
		push(@fe,$1);
	    }
	}
    } else {
	foreach (<OUT>) {
	    if (/^[\(\)\.]+\s+(\-[\d\.]+)/) {
		push(@fe,$1); 
	    }
	}	
    }
    my @stats = (min(\@fe),maxm(\@fe),ave(\@fe),median(\@fe),stdev(\@fe));
    printf "Sample energies: min %.2f, max %.2f, ave %.2f, med %.2f, with stdev %.2f\n",@stats;
    return \@stats;
}

#function to determine if new distribution is plausible wrt original
# currently, if new MFE is within a stdev of old median sample energy
sub plausible {
    my ($stats,$newstats) = @_;
    if ($newstats->[0] < $stats->[3]+$stats->[4]) {
	return 1;
    }
    return 0;
}

sub print_clus {
    my $newclus = shift;
    my $k = 1;
    for (my $i = 1; $i<scalar(@$newclus); $i++) {
	next if (!$newclus->[$i]);
	printf "For cluster %d:\n",$k++;
	my $clus = $newclus->[$i];
	my $feats = $clus->{"constraints"};
	my $coords = $clus->{"seqs"};
	foreach my $seq (keys %$feats) {
	    my $feat = $feats->{$seq};
	    my $coord = $coords->{$seq};
	    print "\t$seq (@$feat) [@$coord]\n";
	}
	my $centroid = $clus->{"centroid"};
	print "  Median centroid: @$centroid\n";
    }
}

sub print_clus_to_file {
    my ($newclus,$fam,$num) = @_;
    my $file = "../../data/analysis/$fam"."_analysis/dbscan_stem.txt";
    $file =~ s/\.txt/$num\.txt/ if $num;
    open CLUS,">$file" or die "cannot open $file\n";
    my $k = 1;
    for (my $i = 1; $i<scalar(@$newclus); $i++) {
	next if (!$newclus->[$i]);
	printf CLUS "For cluster %d:\n",$k++;
	my $clus = $newclus->[$i];
	my $feats = $clus->{"constraints"};
	my $coords = $clus->{"seqs"};
	foreach my $seq (keys %$feats) {
	    my $feat = $feats->{$seq};
	    my $coord = $coords->{$seq};
	    print CLUS "\t$seq (@$feat) [@$coord]\n";
	}
	my $centroid = $clus->{"centroid"};
	print CLUS "  Median centroid: @$centroid\n";
    }
    close(CLUS);
}

#runs metrics on the final clusters to assess 
#precision/recall of consensus stems
sub get_Rfam_ref {
    my ($fam) = @_;
    my $file = "../../data/Rfam_alignmts/$fam". ".txt";
    open RFAM,"<$file" or die "cannot open $file\n";
    my $refseq = "";
    my $native = "";
    foreach (<RFAM>) {
	if (/GC RF\s+([acguACGU\.]+)/) {
	    $refseq = $1;
	} elsif (/SS\_cons\s+([\:\(\)\<\>\[\]\,\.\-\_A]+)/) {
	    $native = $1;
	}
    }
    close(RFAM);
    return [$refseq,$native];
}

#processes rfam reference to get coords of the native bp
sub native_stems {
    my $rfam = shift;
    my @seq = map {$_ eq '.' ? () : ($_)} (split(//,$rfam->[0]));
    my @nat = map {$_ eq '.' ? () : ($_)} (split(//,$rfam->[1]));
    my @bps;
    my $k = 0;
    my $last  = -1;
    my $last2 = -1;
    print "Nat has length ",scalar(@nat),"\n";
    for (my $i = 0; $i < scalar(@nat); $i++) {
	if ($nat[$i] eq '(' || $nat[$i] eq '<' || $nat[$i] eq '[') {
	    push(@bps,$i);
	} elsif ($nat[$i] eq ')' || $nat[$i] eq '>' || $nat[$i] eq ']') {
	    my $left = pop(@bps);
	    if ($last != -1 && ($left != $last - 1 || $last2 != $i-1)) {
		print $last+1," ",$last2+1," $k\n";
		$k = 0;
	    }
	    $last = $left;
	    $last2 = $i;
	    $k++;
	}
    }
    print $last+1," ",$last2+1," $k\n";
}
################################# NOT USED #########################################################

#NOT USED
#makes a dotplot initialized to window 
#all matches get a score of 1, if its part of a HC, it gets its freq added
sub make_dotplot {
    my ($bounds,$nucl) = @_;
    my %dotplot;
    for (my $i=$bounds->[0]; $i <= $bounds->[1]; $i++) {
	my %jdot;
	for (my $j=$bounds->[2]; $j <= $bounds->[3]; $j++) {
	    if (match($nucl->[$i],$nucl->[$j])) {
		$jdot{$j} = 1;
#		print "$nucl->[$i] at $i pairs with $nucl->[$j] at $j\n" if ($i == $bounds->[0]);
	    }
	}
	$dotplot{$i} = \%jdot;
    }
    return \%dotplot;
}

sub print_dotplot {
    my ($dotplot,$bounds) = @_;
    for (my $i=$bounds->[0]; $i <= $bounds->[1]; $i++) {
	if (exists $dotplot->{$i}) {
	    my $jdot = $dotplot->{$i};
	    for (my $j = $bounds->[2]; $j <= $bounds->[3]; $j++) {
		if (exists $jdot->{$j}) { print "$jdot->{$j}  ";}
		else {print "0  ";}
	    }
	} else {
	    my $k = $bounds->[3]-$bounds->[2]+1;
	    my $line = "0  " x $k;
	    print "$line  ";
	}
	print "\n";
    }
}

#calcs a goodness of fit score for each cluster based on conformity to centroid
#calibrate centroid to median length offset
#calculate poisson likelihood of a seq's i,j coordinates:
#calc after i's prob based on needed indels, between i,j's needed indels, and after j's needed indels
sub calc_poisson {
    my ($newclus,$lengths) = @_;
    my $median = $lengths->{"median"};
    my @ratios;
    for (my $i = 1; $i<scalar(@$newclus); $i++) {
	my $clus = $newclus->[$i];
	my $centroid = $clus->{"centroid"};
	my $clusmedian = $clus->{"offset"};
	my $offset = $median - $clusmedian;
	my $medcentroid = [$centroid->[0]+$offset,$centroid->[1]+$offset];
	foreach my $seq (keys %$clus) {
	    next if ($seq eq "centroid" || $seq eq "offset");
	    my $coords = $clus->{$seq};
	    my $length = $lengths->{$seq};
	    my $ratios = calc_single_score($median,$medcentroid,$coords,$length);
	    print "For seq $seq, ratios @$ratios\n";
	    push(@ratios,ave($ratios));
	}
	printf "Average poisson ratio per seq is %.3f\n",ave(\@ratios);
#	last;
    }
}

#calculates the goodness of a seq's hc coordinates based on its expected deviation from median
sub calc_single_poisson {
    my ($median,$medcentroid,$coords,$length) = @_;
    my @ratios;
    my $verbose = 0;
    my $indels = $length - $median; #expected number of indels
#measuring goodness of i coordinate
    my $iseensofar = $coords->[0]-$medcentroid->[0]; #indels seen so far up to i'
    my $left = abs($indels - $iseensofar); #number of indels needed after i'; this is Poisson k
    my $rate = abs($median-$length)/$length; #expected indels per nuc
    my $expected = ceil(($length-$coords->[0])*$rate); #expected indels after i'; this is Poisson lambda
    my $ratio = poisson_ratio($left,$expected);
    push(@ratios,$ratio);
    printf "After  i = $coords->[0], needed $left vs expected $expected: %.3f\n",$ratio if ($verbose);
#measuring goodness of j coordinate in relation to i
    my $jseensofar = $coords->[1]-$medcentroid->[1];
    $left = abs($jseensofar - $iseensofar);
    $expected = ceil(($coords->[1]-$coords->[0])*$rate);
    $ratio = poisson_ratio($left,$expected);
    push(@ratios,$ratio);
    printf "After i $coords->[0] before j $coords->[1], needed $left vs expected $expected: %.3f\n",$ratio if ($verbose);
#measuring goodness of j coordinate
    $left = abs($indels-$jseensofar);
    $expected = ceil(($length-$coords->[1])*$rate);
    $ratio = poisson_ratio($left,$expected);
    push(@ratios,$ratio);
    printf "After j $coords->[1], needed $left vs expected $expected: %.3f\n",$ratio if ($verbose);
    return \@ratios;
}

#calculates the ratio of poisson probabilities
#using p(k,lambda)/p(lambda,lambda)
#equal to L^(k-L)L!/k!
sub poisson_ratio {
    my ($k,$lambda) = @_;
    return 1 if ($k == $lambda);
    my $diff = $k-$lambda;
    my $ratio = 1;
    my $num = $lambda;
    for (my $i = 0; $i < abs($diff); $i++) {
	if ($diff > 0) {
	    if ($k) {
		$ratio *= ($lambda/$k--);
	    } else {
		$ratio *= $lambda;
	    }
	} else {
	    
	    $ratio *= ($num--/$lambda);
	}
    }
    return $ratio;
}
