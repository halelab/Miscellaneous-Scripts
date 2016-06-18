#!/usr/bin/perl

use strict;
use warnings;

print "Enter the name of the input file:\n";
    chomp ( my $in_file = <STDIN> );

unless ( -e $in_file ) {
    print "\n";
    system "dir /B";
    print "\n\nCannot open file! The above files are available for input.\n";
    print "Please re-enter the input file name (with any extensions):\n";
    chomp ( $in_file = <STDIN> );
    unless ( -e $in_file ) {
        print "\nNo such file!\n";
        print "Press enter to exit...\n";
            chomp ( my $exit = <STDIN> );
        exit;
    }
}

open my $in_fh, "<", $in_file;

print "\nWhat is the symbol used for missing data? (e.g. -,.,etc.):\n";
	chomp ( my $missing = <STDIN> );
print "\nDelete markers with more than how many missing data points [RILs]? (e.g. 1,5,10,etc.):\n";
	chomp ( my $missing_threshold = <STDIN> );
print "\nDelete markers if and only if an acceptable marker is within how many cM? (e.g. 0.5,1,1.5,etc.):\n";
	chomp ( my $missing_neighbor = <STDIN> );
print "\nSearch for DCO's within a window length of what size [cM]? (e.g. 0.5,1,1.5,etc.):\n";
	chomp ( my $window_length = <STDIN> );
print "\nChoose the size [cM] of the scanning window (e.g. 1,2,5,etc.):\n";
	chomp ( my $scan_window = <STDIN> );

print "\nProvide a filename for the tab-delimited results:\n";
	chomp ( my $results_file = <STDIN> );

if ( -e $results_file && -s _ ) {
    print "\n'$results_file' already exists, overwrite? [y/n]\n";
    chomp ( my $overwrite = <STDIN> );
    if ( $overwrite =~ /[nN]/ ) {
        print "\nProvide another filename for the results:\n";
        chomp ( $results_file = <STDIN> );
        if ( -e $results_file && -s _ ) {
            print "\n'$results_file' already exists. Press enter to exit.\n";
            chomp ( my $exit = <STDIN> );
            exit;
        }
    }
}

my @map = ();

while ( <$in_fh> ) {
	chomp;
	push @map, [ split /,/ ];
}

close $in_fh;

my $row1 = $map[0];
my $total_RILs = scalar ( @$row1 ) - 5;
my $total_markers = scalar ( @map );
my $removed_due_to_missing = 0;
my $removed_due_to_DCO = 0;

my @reduced_map = ();
push @reduced_map, [ @$row1 ];

# Remove markers based on missing data (excluding the most distal markers on each chromosome)

for ( my $j=1; $j <= $total_markers - 2; $j++ ) {
	my $ref_row = $map[$j];
	
	if ( $map[$j-1][1] ne $map[$j][1] || $map[$j][1] ne $map[$j+1][1] ) {
		push @reduced_map, [ @$ref_row ];
		next;
	}
	
	if ( $map[$j][2] > $missing_threshold &&
			( ( ( $map[$j+1][4] - $map[$j][4] ) < $missing_neighbor && $map[$j+1][2] <= $missing_threshold ) ||
			  ( ( $map[$j][4] - $map[$j-1][4] ) < $missing_neighbor && $map[$j-1][2] <= $missing_threshold ) ) ) {
		$removed_due_to_missing++;
		next;
	} else {
		push @reduced_map, [ @$ref_row ];
		next;
	}

}

my $row_last= $map[ $total_markers - 1 ];
push @reduced_map, [ @$row_last ];

my $removed_markers = $removed_due_to_missing;

print "\nMarkers removed due to missing data = $removed_markers\n";

# Remove markers causing simple DCO's (ABA or BAB) in the map

print "\nBegin iterative rounds to remove simple DCO's (window length = $window_length):\n\n";

my $meta_convergence = 0;
my $convergence = 0;
my $frame1_convergence = 0;
my $frame2_convergence = 0;
my $frame3_convergence = 0;
my @reduced2_map = ();

while ( $meta_convergence == 0 ) {

	while ( $convergence == 0 ) {

		$removed_due_to_DCO = 0;

		@reduced2_map = ();
		my $total_reduced_markers = scalar ( @reduced_map );

		for ( my $j=1; $j <= $total_reduced_markers - 2; $j=$j+3 ) {

			my $pre_row = $reduced_map[$j-1];
			my $ref_row = $reduced_map[$j];
			my $post_row = $reduced_map[$j+1];
		
			push @reduced2_map, [ @$pre_row ];

			if ( @$pre_row[1] ne @$ref_row[1] || @$ref_row[1] ne @$post_row[1] ) {
				push @reduced2_map, [ @$ref_row ];
				push @reduced2_map, [ @$post_row ];
				next;
			}

			my $DCO_count = 0;
			for ( my $i=5; $i <= $total_RILs + 4; $i++ ) {
				if ( 	$reduced_map[$j-1][$i] eq $reduced_map[$j+1][$i] &&
						$reduced_map[$j-1][$i] ne $reduced_map[$j][$i] &&
						$reduced_map[$j-1][$i] ne $missing &&
						$reduced_map[$j][$i] ne $missing &&
						$map[$j+1][4] - $map[$j-1][4] <= $window_length ) {
					$DCO_count++;
					next;
				} else {
					next;
				}
			}

			if ( $DCO_count == 0 ) {
				push @reduced2_map, [ @$ref_row ];
				push @reduced2_map, [ @$post_row ];
				next;
			} else {
				$removed_due_to_DCO++;
				push @reduced2_map, [ @$post_row ];
				next;
			}

		}

		my $frame_offset = ($total_reduced_markers) % 3;
	
		if ( $frame_offset == 1 ) {
			my $last_row = $reduced_map[ $total_reduced_markers - 1];
			push @reduced2_map, [ @$last_row ];
		}
	
		if ( $frame_offset == 2 ) {
			my $penult_row = $reduced_map[ $total_reduced_markers - 2];
			push @reduced2_map, [ @$penult_row ];
			my $last_row = $reduced_map[ $total_reduced_markers - 1];
			push @reduced2_map, [ @$last_row ];
		}
	
		@reduced_map = @reduced2_map;
		$removed_markers = $removed_markers + $removed_due_to_DCO;

		print "Markers removed due to simple DCO's = $removed_due_to_DCO\n";

		if ( $removed_due_to_DCO == 0 ) {
			$convergence = 1;
		} else {
			$convergence = 0;
			$frame1_convergence++;
		}

	}

	# Shift one base and keep iterating...

	$convergence = 0;

	while ( $convergence == 0 ) {

		$removed_due_to_DCO = 0;

		@reduced2_map = ();
		my $total_reduced_markers = scalar ( @reduced_map );

		my $first_row = $reduced_map[0];
		push @reduced2_map, [ @$first_row ];

		for ( my $j=2; $j <= $total_reduced_markers - 2; $j=$j+3 ) {

			my $pre_row = $reduced_map[$j-1];
			my $ref_row = $reduced_map[$j];
			my $post_row = $reduced_map[$j+1];
		
			push @reduced2_map, [ @$pre_row ];

			if ( @$pre_row[1] ne @$ref_row[1] || @$ref_row[1] ne @$post_row[1] ) {
				push @reduced2_map, [ @$ref_row ];
				push @reduced2_map, [ @$post_row ];
				next;
			}

			my $DCO_count = 0;
			for ( my $i=5; $i <= $total_RILs + 4; $i++ ) {
				if ( 	$reduced_map[$j-1][$i] eq $reduced_map[$j+1][$i] &&
						$reduced_map[$j-1][$i] ne $reduced_map[$j][$i] &&
						$reduced_map[$j-1][$i] ne $missing &&
						$reduced_map[$j][$i] ne $missing &&
						$map[$j+1][4] - $map[$j-1][4] <= $window_length ) {
					$DCO_count++;
					next;
				} else {
					next;
				}
			}

			if ( $DCO_count == 0 ) {
				push @reduced2_map, [ @$ref_row ];
				push @reduced2_map, [ @$post_row ];
				next;
			} else {
				$removed_due_to_DCO++;
				push @reduced2_map, [ @$post_row ];
				next;
			}

		}

		my $frame_offset = ($total_reduced_markers - 1) % 3;
	
		if ( $frame_offset == 1 ) {
			my $last_row = $reduced_map[ $total_reduced_markers - 1];
			push @reduced2_map, [ @$last_row ];
		}
	
		if ( $frame_offset == 2 ) {
			my $penult_row = $reduced_map[ $total_reduced_markers - 2];
			push @reduced2_map, [ @$penult_row ];
			my $last_row = $reduced_map[ $total_reduced_markers - 1];
			push @reduced2_map, [ @$last_row ];
		}
	
		@reduced_map = @reduced2_map;
		$removed_markers = $removed_markers + $removed_due_to_DCO;

		print "Markers removed due to simple DCO's = $removed_due_to_DCO\n";

		if ( $removed_due_to_DCO == 0 ) {
			$convergence = 1;
		} else {
			$frame2_convergence++;
			$convergence = 0;
		}

	}

	# Shift one base and keep iterating...

	$convergence = 0;

	while ( $convergence == 0 ) {

		$removed_due_to_DCO = 0;

		@reduced2_map = ();
		my $total_reduced_markers = scalar ( @reduced_map );

		my $first_row = $reduced_map[0];
		my $second_row = $reduced_map[1];
		push @reduced2_map, [ @$first_row ];
		push @reduced2_map, [ @$second_row ];

		for ( my $j=3; $j <= $total_reduced_markers - 2; $j=$j+3 ) {

			my $pre_row = $reduced_map[$j-1];
			my $ref_row = $reduced_map[$j];
			my $post_row = $reduced_map[$j+1];
		
			push @reduced2_map, [ @$pre_row ];

			if ( @$pre_row[1] ne @$ref_row[1] || @$ref_row[1] ne @$post_row[1] ) {
				push @reduced2_map, [ @$ref_row ];
				push @reduced2_map, [ @$post_row ];
				next;
			}

			my $DCO_count = 0;
			for ( my $i=5; $i <= $total_RILs + 4; $i++ ) {
				if ( 	$reduced_map[$j-1][$i] eq $reduced_map[$j+1][$i] &&
						$reduced_map[$j-1][$i] ne $reduced_map[$j][$i] &&
						$reduced_map[$j-1][$i] ne $missing &&
						$reduced_map[$j][$i] ne $missing &&
						$map[$j+1][4] - $map[$j-1][4] <= $window_length ) {
					$DCO_count++;
					next;
				} else {
					next;
				}
			}

			if ( $DCO_count == 0 ) {
				push @reduced2_map, [ @$ref_row ];
				push @reduced2_map, [ @$post_row ];
				next;
			} else {
				$removed_due_to_DCO++;
				push @reduced2_map, [ @$post_row ];
				next;
			}

		}

		my $frame_offset = ($total_reduced_markers - 2) % 3;
	
		if ( $frame_offset == 1 ) {
			my $last_row = $reduced_map[ $total_reduced_markers - 1];
			push @reduced2_map, [ @$last_row ];
		}
	
		if ( $frame_offset == 2 ) {
			my $penult_row = $reduced_map[ $total_reduced_markers - 2];
			push @reduced2_map, [ @$penult_row ];
			my $last_row = $reduced_map[ $total_reduced_markers - 1];
			push @reduced2_map, [ @$last_row ];
		}
	
		@reduced_map = @reduced2_map;
		$removed_markers = $removed_markers + $removed_due_to_DCO;

		print "Markers removed due to simple DCO's = $removed_due_to_DCO\n";

		if ( $removed_due_to_DCO == 0 ) {
			$convergence = 1;
		} else {
			$frame3_convergence++;
			$convergence = 0;
		}

	}

	if ( $frame1_convergence == 0 && $frame2_convergence == 0 && $frame3_convergence == 0 ) {
		$meta_convergence = 1;
	} else {
		$meta_convergence = 0;
		$frame1_convergence = 0;
		$frame2_convergence = 0;
		$frame3_convergence = 0;
	}

}

# Output summary stats and write to output file

my $final_marker_count = scalar (@reduced_map);
my $total_DCO_count = $total_markers - $final_marker_count - $removed_due_to_missing;
my $total_removed = $total_markers - $final_marker_count;

print "\nMarkers removed due to DCO's = $total_DCO_count\n";

print "\nSTAGE 1 summary:\n";
print "Total initial markers = $total_markers\n";
print "Markers removed due to missing data and DCO's = $total_removed\n";
print "Remaining marker count = $final_marker_count\n\n";
print "Using the specified window size to select representative markers";

# Select "best" representative marker within each window

my @master = @reduced_map;
my @out_map = @reduced_map;

my $master_markers = scalar ( @master );
my %final_list;

while ( ( scalar ( @out_map ) ) > 0 ) {

	my @map = @out_map;
	my %group_linked;
	my %group_missing;
	my $chromosome = $map[0][1];
	my $start_pos = $map[0][3];
	
	print ".";
	
	for ( my $j=0; $j < scalar ( @map ); $j++ ) {
		if ( ( $map[$j][3] - $start_pos ) <= $scan_window && $map[$j][1] eq $chromosome ) {
			$group_linked{$map[$j][0]} = $map[$j][4];
			$group_missing{$map[$j][0]} = $map[$j][2];
			next;
		} else {
			last;
		}
	}
		
	my @rank = ();
		
	foreach my $locus ( sort { $group_linked{$b} <=> $group_linked{$a} } keys %group_linked ) {
		push @rank, $locus;
	}
	
	my %reduced_group;
	my $linked_tally = 0;
	
	for ( my $j=0; $j < scalar ( @rank ); $j++ ) {
		if ( $group_linked{$rank[$j]} >= $linked_tally ) {
			$linked_tally = $group_linked{$rank[$j]};
			$reduced_group{$rank[$j]} = $group_missing{$rank[$j]};
			next;
		} else {
			last;
		}
	}
		
	my @missing_rank = ();
	
	foreach my $locus ( sort { $reduced_group{$a} <=> $reduced_group{$b} } keys %reduced_group ) {
		push @missing_rank, $locus;
	}
	
	$final_list{$missing_rank[0]} = $missing_rank[0];
		
	@out_map = ();
	
	for ( my $j=0; $j < scalar ( @map ); $j++ ) {
		my $row = $map[$j];
		if ( $group_linked{$map[$j][0]} ) {
			next;
		} else {
			push @out_map, [ @$row ];
			next;
		}
	}
}

# Write to output file

my $end_marker_count = 0;

open my $results_fh, ">", $results_file;

for ( my $j=0; $j <= $master_markers - 1; $j++ ) {
	my $row = $master[$j];
	if ( $final_list{$master[$j][0]} ) {
		my $line = join(',', @$row);
	    print $results_fh "$line\n";
		$end_marker_count++;
		next;
	} else {
		next;
	}
}

close $results_fh;

print "\n\nSTAGE 2 summary:\n";
print "Final marker count = $end_marker_count\n\n";

exit;
