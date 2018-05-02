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
    manpage();
} else {
    my $args = process_args(@ARGV);
    my $verbose = $args->{"verbose"};
    my $type = "Sfold";
    my $names = run_sfold_profiling($args);
#if you've already run Sfold and don't want to rerun every time, comment above and uncomment next two lines
#    my @names = map {get_name($_)} @{$args->{"seqs"}};
#    my $names = \@names;
    my $newclus = 0;
    my $feats;
    my $repeat = 1;
    my $found = {};
    my $extend = {};
    while ($repeat) {
	my $data = process_profiling($args->{"output"},$names,$extend);
	$feats = $data->[0];
	my $lengths = $data->[1];
	my $cells = process_trip_feat($feats,$lengths,$found,$verbose);
	my $eps = find_eps($lengths,$verbose);
	my $clusters = cluster($cells,$eps,$verbose);
	my $clusstems = print_clusters($clusters->[0],$cells,$feats,$verbose);
	$newclus = find_centroid($clusstems->[0],$lengths,$data->[2],$verbose);
#	print_clus_stats($newclus);
	improve_cluster($newclus,$data,$args,$clusstems->[1],$clusters->[1]);
	if (found_newclus($found,$newclus,$feats,$verbose)) {
	    $found = get_truepos($newclus,$feats);
#	    print_clus_stats($newclus);
	    $extend = resample($found,$args,$feats);
	    print "Resampling:\n" if ($verbose);
	} else {
	    $repeat = 0;
	}
    }
    print "FINAL CLUSTERS:\n";
#    print_clus_stats($newclus);
    print_clus($newclus);
    print_clus_to_file($newclus,$args);
}

#usage instructions 
sub manpage {
    print "ConsensusStems code for Academic Users, Version 1.0 (April 2018)\n\n";
    print "Usage: perl ConsensusStems.pl [options]... directory_of_sequences\n\n";
    print "Required: perl, Boltzmann sampling with constraints, RNA profiling, a MFE prediction program\n";
    print "\tBoltzmann sampling with constraints is implemented by Sfold, which can be requested from http://sfold.wadsworth.org/cgi-bin/index.pl\n";
    print "\tRNA profiling is available at https://github.com/gtfold/RNAStructProfiling\n";
    print "\tA recommended MFE prediction program is GTfold's gtmfe option, available at https://github.com/gtfold/gtfold\n\n";
    print "Input a directory containing all sequences (in FASTA format) to be processed.\n\tRecommended: at least 8 related sequences.\n\tConsensusStems is not designed for nor tested with less than 8 sequences\n\n";
    print "Final output with clusters written to [output directory]/ConsensusStems_output.txt\n";
    print "Options:\n";
    print "\t-s <path>\t Path to invoke Sfold, eg -s /usr/bin/sfold (default)\n";
    print "\t-p <path>\t Path to invoke RNAprofiling, eg -p ../RNAprofiling (default is ./RNAprofile)\n";
    print "\t-m <path>\t Path to invoke MFE program, eg -m /usr/bin/gtmfe (default)\n";
    print "\t-d <name>\t Path to MFE NNTM parameters eg -d Desktop/data/Turner99 (default is ./data)\n";
    print "\t-o <path>\t Path to store output files, eg -o ../output (default is ./output)\n";
    print "\t-f <name>\t Family name, eg -f tRNA (default is last directory name if exists or RNA)\n";
    print "\t-v \t Flag to enable verbose output concerning details of the method\n";

}

sub process_args {
    my @args = @_;
    my $verbose = 0;
    my $arghash = initialize();
    my $dir = pop(@args);
    $arghash->{"seqs"} = process_seqs($dir);
    while (scalar(@args) > 0) {
	process_options($arghash,\@args);
    }
    if (!exists $arghash->{"family"}) {
	my @dirs = split(/\//,$dir);
	if (scalar(@dirs)>1) {
	    $arghash->{"family"} = pop(@dirs);
	} else {
	    $arghash->{"family"} = "RNA";
	}
    }
    if ($verbose) {
	print "Sfold path is ",$arghash->{"sfold"},"\n";
	print "Prof path is ",$arghash->{"prof"},"\n";
	print "MFE path is ",$arghash->{"mfe"},"\n";
	print "output path is ",$arghash->{"output"},"\n";
	print "family is ",$arghash->{"family"},"\n";
    }
    return $arghash;
}

sub initialize {
    my %args;
    $args{"sfold"} = "/usr/bin/sfold";
    $args{"mfe"} = "/usr/bin/gtmfe";
    $args{"prof"} = "./RNAprofile";
#    $args{"family"} = "RNA";
    $args{"output"} = "./output/";
    $args{"verbose"} = 0;
    return \%args;
}

sub process_seqs {
    my $dir = shift;
    my @dir = split(//,$dir);
#    print "Found dir ",@dir," with last char $dir[-1]\n";
    $dir = "$dir"."/" if ($dir[-1] ne '/');
    my $files = `ls $dir | wc -l`;
    print "Warning! Found less than eight files!\n" if ($files < 8);
    my @files = `ls $dir`;
    my @newfiles;
    foreach my $file (@files) {#one file per line assumed
	chomp($file);
	$file = "$dir" . "$file";
	open SEQ, "$file" or die "cannot open $file\n";
	my $first = 1;
	foreach (<SEQ>) {
	    die "File $file is not in FASTA format\n" if ($first && !/^\>/);
	    last;
	}
	push(@newfiles,$file);
    }
    return \@newfiles;
}

sub process_options {
    my ($arghash,$args) = @_;
    while ($args->[0] ne "-s" && $args->[0] ne "-p" && $args->[0] ne "-f" && $args->[0] ne "-o" && $args->[0] ne "-m" && $args->[0] ne "-d" && $args->[0] ne "-v") {
	shift(@$args);
	return if (scalar(@$args)<1);
    }
    if ($args->[0] eq "-s") {
	$arghash->{"sfold"} = $args->[1];
	splice(@$args,0,2);
    } elsif ($args->[0] eq "-p") {
	$arghash->{"prof"} = $args->[1];
	splice(@$args,0,2);
    } elsif ($args->[0] eq "-m") {
	$arghash->{"mfe"} = $args->[1];
	splice(@$args,0,2);
    } elsif ($args->[0] eq "-d") {
	$arghash->{"paramdir"} = $args->[1];
	splice(@$args,0,2);
    } elsif ($args->[0] eq "-f") {
	$arghash->{"family"} = $args->[1];
	splice(@$args,0,2);
    } elsif ($args->[0] eq "-v") {
	$arghash->{"verbose"} = 1;
	shift(@$args);
    } elsif ($args->[0] eq "-o") {
	my @out = split(//,$args->[1]);
	$args->[1] = "$args->[1]" . "/" if ($out[-1] ne '/');
	$arghash->{"output"} = $args->[1];
	splice(@$args,0,2);
    }
}


#runs profiling on sfold outputs of fam
sub run_sfold_profiling {
    my ($args) = @_;
    my $fam = $args->{"family"};
    my $seqs = $args->{"seqs"};
    my $spath = $args->{"sfold"};
    my $ppath = $args->{"prof"};
    my $output = $args->{"output"};
    `mkdir $output`;
    my @names=();
    foreach my $seqfile (@$seqs) {
	my $name = get_name($seqfile);
	push(@names,$name);
	my $profile = "$output"."$name.out";
	my $outdir = "$output"."Sfold"."_$name";
	print "Running Sfold on $seqfile\n";
	`$spath -o $outdir $seqfile`;
	my $outfile = $outdir."/sample_1000.out";
	`$ppath -sfold $outfile -g -v $seqfile > $profile`;
#	last;
    }
#    $args->{"seqs"} = \@names;
    return \@names;
}

sub get_name {
    my $seq = shift;
    my @dirs = split(/\//,$seq);
    $seq = $dirs[-1];
    $seq = $1 if ($seq =~ /(\S+)\.txt/);
    $seq = $1 if ($seq =~ /(\S+)\.fa/);
    $seq = $1 if ($seq =~ /(\S+)\.fasta/);
    $seq = $1 if ($seq =~ /(\S+)\.seq/);
    return $seq;
}

#grabs feature triplet info and puts into hash
#processes input profiling file based on type:
# 0 = gtboltzmann input, 1 = sfold input, 2 = resample input
sub process_profiling {
    my ($outdir,$names,$extension) = @_;
    my %feats;
    my %lengths;
    my %nucl;
    my %files;
    foreach my $seq (@$names) {
	my $extend = ".out";
	if (exists $extension->{$seq}) {$extend = "_resample.out";}
	my $file = "$outdir"."$seq"."$extend";
	my @feat = ();
	open OUT, "<$file" or die "cannot open $file\n";
	foreach (<OUT>) {
	    if (/Featured helix (\d+)\: (\d+ \d+ \d+) with freq (\d+)/) {
#	    $feats{$1} = $2;
		$feat[$1] = $2;
	    } elsif (/seq in .+ is ([UCAGTucagt]+) with length (\d+)/) {
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


#compare against already found features in %found; add if not present
sub process_trip_feat {
    my ($feats,$lengths,$found,$verbose) = @_;
    my $median = median([values %$lengths]);
    my %featcells;
#    my @heatmap = (0) x $median;
    foreach my $seq (keys %$feats) {
	print "processing $seq\n" if ($verbose);
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
		print "Adding $seq ($k) $coord\n" if ($verbose);
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
    my ($lengths,$verbose) = @_;
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
	print "found eps $eps with $noise seq as noise\n" if ($verbose);
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
    my ($cells,$params,$verbose) = @_;
    my $eps = $params->[0];
    my $minpts = $params->[1];
    my $n = scalar(keys %$cells);
    my $noise = $n;
    my $label;
    while ($noise > $n/2) {
	if ($minpts < 3) {
	    print "minpts below critical level of 3; raising eps to ",$eps+1,"\n" if ($verbose);
	    $minpts = $params->[1];
	    $eps = 3;
	}
	$label = DBscan($cells,$eps,$minpts,"coord");
	my @noise = map {$label->{$_} == 0 ? ($_):() } keys %$label;
	my @n = map {scalar(@{$cells->{$_}}) } @noise;
	$noise = sum(\@n);
	print "found minpts $minpts with $noise points as noise out of $n\n" if ($verbose);
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
    my ($label,$cells,$stems,$verbose) = @_;
    my %clus;
    my %found;
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
	}
	my @sorted = sort {$a cmp $b} keys %allstems;
	my %myclus;
	foreach my $stem (@sorted) { #for every stem in the cluster
	    my @stem = split(/\s/,$stem);
	    my $steminfo;
	    my $steminfoprint;
	    my $coords = $stems->{$stem[0]}[$stem[1]];
	    $steminfoprint = "($stem[1]) [$coords]";
	    $steminfo = "$stem[1] | $coords";
	    print "\t$stem[0] $steminfoprint\n" if ($verbose);
	    add(\%myclus,$stem[0],$steminfo);
	    $found{$stem} = $clust if ($clust);
	}
#	map {$found{$_} = $clust} @sorted if ($clust);
	$clusters[$clust] = \%myclus;
    }
    return [\@clusters,\%found];
}


#combines all stems of a seq into a superstem
#finds the median coordinates of all seq stems
#variation: use normalized coords?
sub find_centroid {
    my ($clusters,$lengths,$nucls,$verbose) = @_;
    my @newclus;
    my %hctoclus;
    my $median = $lengths->{"median"};
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
	    my $bounds =  [$i,$i+$k-1,$j-$l+1,$j];
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
    my ($newclus,$data,$args,$found,$minpts) = @_;
    my $feats = $data->[0];
    my $lengths = $data->[1];
    my $median = $lengths->{"median"};
    my @seqs = keys %$lengths;
    remove_char(\@seqs,"median");
    my %combine;
    my $verbose = $args->{"verbose"};
#make dir to store mfe seqs for missing seqs
    my $outdir = $args->{"output"};
    my $condir = $outdir . "constraints";
    `mkdir $condir`;
    $outdir = $outdir . "seqs";
    `mkdir $outdir`;
#processing each cluster
    for (my $i = 1; $i<scalar(@$newclus); $i++) {
	my $clus = $newclus->[$i];
	my $size = scalar(keys %{$clus->{"seqs"}});
	if ($size < $minpts) {
#	    my $badclus = $newclus->[$i];
	    $newclus->[$i] = 0;
	    print "\tCluster $i with ",$size," below $minpts threshold: deleting\n" if ($verbose);
	    next;
	}
#make dir for each cluster
	my $clusdir = $outdir."/cluster_$i";
	`mkdir $clusdir`;
	my ($score,$newconstraints) = @{find_missing(@_,$i)};
	my $total = $score + $size; #since each present member of clus gets score of 1, so sum = size
	if ($total <= 0) {
	    $newclus->[$i] = 0;
	    print "\tCluster $i with score of $total: deleting\n" if ($verbose);
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
		print "\tAdding $seq (@{$constraints->{$seq}}) [@{$stemcoords->{$seq}}]\n" if ($verbose);
	    }
	    my $newcentroid = recalc_centroid($stemcoords);
	    $clus->{"centroid"} = $newcentroid;
	    print "\tNew centroid: @$newcentroid\n" if ($verbose);
	    print "\tTotal score of cluster $i: $total\n" if ($verbose);
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
    my ($newclus,$data,$args,$found,$minpts,$i) = @_;
    my $feats = $data->[0];
    my $lengths = $data->[1];
    my $clus = $newclus->[$i];
    my $median = $lengths->{"median"};
    my $missing = $clus->{"missing"};
    my $centroid = $clus->{"centroid"};
    my $bpnums = $clus->{"bpnum"};
    my $med_coords = $clus->{"medcoords"};
    my $verbose = $args->{"verbose"};
    my $totalscore = 0;
    my %constrain;
    print "For cluster $i:\n" if ($verbose);
    foreach my $seq (@$missing) {
	print "   For seq $seq, " if ($verbose);
	my $feat = $feats->{$seq};
	my $length = $lengths->{$seq};
	my $diff = $length-$median;
	my $nucl = $data->[2]{$seq};
	my $bounds = find_boundaries($med_coords,$diff,$length);
	my $scores = search_window($seq,$bounds,$nucl,$centroid,$found,$diff,$args,$newclus,$i,$lengths,$data->[3]);
	my $mfe = run_mfe($args,$i,$seq,$nucl,$bounds);
	my $bpnum = 0;
	if ($mfe->[1] != 0) {
	    foreach my $bp (@{$mfe->[0]}) {
		if ($bp eq '(') {$bpnum++;}
	    }
	    if ($bpnum >= $bpnums) {$totalscore++;}
#	    elsif ($bpnum >= ceil($bpnums/2)) {$scores{$seq} = 0;}
	    elsif ($bpnum < ceil($bpnums/2)) {$totalscore--;}
	} else {
	    $totalscore -= 2;
	}

#if the found is a feature, save it
	my @found = map {$_->[0] < scalar(@$feat) ? ($_->[0]) : ()} @$scores;
	if (scalar(@found)<1) {next;}
	$constrain{$seq} = \@found;
	print "found potential features @found\n" if ($verbose);
    }
    print "Cumulative score of missing seqs: $totalscore\n" if ($verbose);
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
    my ($args,$i,$seq,$nucl,$bounds) = @_;
    my $outdir = $args->{"output"};
    my $fam = $args->{"family"};
    my $verbose = 0;
    my $mkdirs = 1;
#    my $name = "$outdir" . "$seq";    
#make sequence file to run MFE
    my $file = "$outdir"."seqs/cluster_$i"."/$seq".".txt";
    open SEQ, ">$file" or die "cannot open $file\n";
    print SEQ ">$fam, $seq for window @$bounds\n";
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
    my $constraints = $outdir."constraints/$seq" . ".txt";
    make_constraints($constraints,$diff,$k,$length);
    my $gtmfe = $args->{"mfe"};
    my $paramdir = $args->{"paramdir"};
    my @mfe = `$gtmfe --paramdir $paramdir -c $constraints $file`;
    print "@mfe\n" if $verbose;
    my $mfe = process_mfe($diff,$k,@mfe);
    return $mfe;
}

#write different code to process mfe output here
#if gtmfe format not used
sub process_mfe {
    my ($diff,$k,@mfe) = @_;
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
    return [$bp,$mfe];
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
    my ($seq,$bounds,$nucl,$centroid,$found,$diff,$fam,$allclus,$i,$lengths,$files) = @_;
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
    my ($found,$newclus,$feats,$verbose) = @_;
    return 1 if (scalar(keys %$found) == 0);
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
    my ($found,$args,$feats) = @_;
    my $seqs = $args->{"seqs"};
    my $out = $args->{"output"};
    my $spath = $args->{"sfold"};
    my $ppath = $args->{"prof"};
    my $verbose = $args->{"verbose"};
    my %found;
    my $mkdir = 1;
    my %extend;
#prohibit all features not found in clusters
    foreach my $seqfile (@$seqs) {
	my $seq = get_name($seqfile);
	my $feat = $feats->{$seq};
	my $myfound = [];
	$myfound = $found->{$seq} if (exists $found->{$seq});
	print "examining $seq with found @$myfound\n" if ($verbose);
	if (scalar(@$feat) <= scalar(@$myfound)+1) {next;}
	my $cfile = "$out"."$seq"."_constraints.txt";
	open CON, ">$cfile" or die "cannot open $cfile\n";	
	for (my $i = 1; $i<scalar(@$feat); $i++) {
	    if (exists_char($myfound,$feat->[$i])) { next;}
	    print CON "P $feat->[$i]\n";
	}
	close(CON);
#run sfold with constraints
	my $outdir = "$out"."$seq"."_resample";
	print "Resampling $seqfile\n";
	`$spath -o $outdir -f $cfile $seqfile`;
#test new structures are energetically plausible
	my $sout = "$out"."Sfold"."_$seq"."/sample_1000.out";
	my $stats = energy_calcs($sout);
	$sout = $outdir . "/sample_1000.out";
	my $newstats = energy_calcs($sout,$verbose);
	if (!plausible($stats,$newstats)) {
	    print "implausible new sample for $seq; sticking with original\n" if ($verbose);
	    next;
	}
#run profiling on new sfold output
	my $profile = "$out"."$seq"."_resample.out";
	`$ppath -sfold $sout -g -v $seqfile > $profile`;
	$extend{$seq}++;
    }
    return \%extend;
}

sub energy_calcs {
    my ($out,$verbose) = @_;
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
    printf "Sample energies: min %.2f, max %.2f, ave %.2f, med %.2f, with stdev %.2f\n",@stats if ($verbose);
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


sub print_clus_stats {
    my $newclus = shift;
    for (my $i = 1; $i<scalar(@$newclus); $i++) {
	next if (!$newclus->[$i]);
	my $clus = $newclus->[$i];
	my $feats = $clus->{"constraints"};
	printf "found %d seqs for cluster %d\n", scalar(keys %$feats),$i;
    }
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
    my ($newclus,$args) = @_;
    my $outdir = $args->{"output"};
    my $file = $outdir . "ConsensusStems_output.txt";
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
