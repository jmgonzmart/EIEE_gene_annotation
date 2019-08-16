#!/software/bin/perl -w

#This script compares transcript annotation in a gene list between different releases

#Counts of:
#-novel exons
#-novel introns
#-shifted splice junctions
#by gene and in total

#-genomic coverage by gene and for all new features (needs to be additive i.e. don’t subtract bases from the count for a shifted splice site that creates a shorter exon)

#Print tsv as follows - one line per transcript - include transcripts with no reported novel features:
#-Gene symbol
#-Ensembl gene ID
#-Ensembl transcript ID
#-Otter gene ID
#-Otter transcript ID
#-transcript biotype
#-number of novel exons
#-number of novel introns
#-number of shifted splice junctions

#And finally at the level of gene models - GTF/GFF for novel features
#-exon trios plus flanking introns for novel exons
#-exon pairs flanking intron for novel introns
#-exon pairs flanking intron for shifted splice junctions


#Exclude retained_intron transcripts except for transcript count


use strict;
use Getopt::Long;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Feature;

$|=1;

my $infile;
my $outfile;
my $species;

&GetOptions(
            'in=s'      => \$infile,
            'out=s'     => \$outfile,
            'species=s' => \$species,
           );



#Old annotation
my $old_annot_db; 
if ($species eq "human"){
  $old_annot_db = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
    -user   => 'anonymous',
    -dbname => 'homo_sapiens_core_76_38', #Gencode 20
    -host   => 'ensembldb.ensembl.org',
    -driver => 'mysql'
  );
}
elsif ($species eq "mouse"){
  $old_annot_db = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
    -user   => 'anonymous',
    -dbname => 'mus_musculus_core_74_38', #Gencode  M2
    -host   => 'ensembldb.ensembl.org',
    -driver => 'mysql'
  );
}
my $oa_sa = $old_annot_db->get_SliceAdaptor;
my $oa_ga = $old_annot_db->get_GeneAdaptor;
my $oa_ta = $old_annot_db->get_TranscriptAdaptor;


#New annotation
my $new_annot_db;
my $core_db; #To get ENS tr ids when loutre is used
if ($species eq "human"){
  $new_annot_db = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
    -user   => 'anonymous',
    -dbname => 'homo_sapiens_core_92_38', #Gencode 28
    -host   => 'ensembldb.ensembl.org',
    -driver => 'mysql'
  );
}
elsif ($species eq "mouse"){
  $new_annot_db = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
    -user   => 'anonymous',
    -dbname => 'mus_musculus_core_90_38', #Gencode M15
    -host   => 'ensembldb.ensembl.org',
    -driver => 'mysql'
  );
}
my $na_sa = $new_annot_db->get_SliceAdaptor;
my $na_ga = $new_annot_db->get_GeneAdaptor;
my $na_ta = $new_annot_db->get_TranscriptAdaptor;
my $c_ta = $core_db->get_TranscriptAdaptor if $core_db;

open (OUT, ">$outfile") or die "Can't open $outfile\n";
print OUT "GENE\tGene symbol\tEnsembl gene id\tNew exons\tNew introns\tExon skipping events\tShifted SJs\tNew exon coverage\tShifted SJ coverage\tExtra CDS coverage\tExtra exon bases\tExtra CDS bases\tExtra 5' UTR bases\tExtra 3' UTR bases\n";
print OUT "TR\tGene symbol\tEnsembl gene id\tEnsembl transcript id\tTranscript biotype\tStart/end_NF tags\tNew exon type\tCreated date\tNew exons\tNew introns\tExon skipping events\tShifted SJs\n";

open (IDS4GTF, ">$outfile.ids4gtf") or die "Can't open $outfile.ids4gtf\n";
my $c = 0;
open (LIST, $infile) or die "Can't open $infile\n";
while (<LIST>){
  chomp;
  my ($ensid, $gsymbol) = split(/\t/);
  my $gene_o = $oa_ga->fetch_by_stable_id($ensid);# or die "Can't find old gene $ensid";
  if (!$gene_o){
    foreach my $gene (@{$oa_ga->fetch_all_by_external_name($gsymbol)}){
      if ($gene->is_reference){
        $gene_o = $gene;
        print "Using ".$gene_o->stable_id." as old gene for $gsymbol\n";
      }
    }
  }
  die "Can't find old gene $ensid / $gsymbol" unless $gene_o;
  my $gene_n = $na_ga->fetch_by_stable_id($ensid);# or die "Can't find new gene $ensid";
  if (!$gene_n){
    foreach my $gene (@{$na_ga->fetch_all_by_external_name($gsymbol)}){
      if ($gene->is_reference){
        $gene_n = $gene;
        print "Using ".$gene_n->stable_id." as new gene for $gsymbol\n";
      }
    }
  }
  die "Can't find new gene $ensid / $gsymbol" unless $gene_n;
  
  print OUT "\n\n#############\n\nChecking gene $ensid / $gsymbol\n";



#two tdf files with lists of
#1 - all v20 transcripts
#2 - all transcripts in new but not v20
#And columns:
# Gene symbol
# Ensembl gene ID
# Ensembl transcript ID


#   #by OTT id or by exon structure?
#   my %old_trids;
#   my %old_exon_str;
#   foreach my $tr_o (@{$gene_o->get_all_Transcripts}){ 
#     my $old_ott_tid = Gencode::Default->get_otter_id($tr_o);
#     $old_trids{$old_ott_tid}++ if $old_ott_tid;
#     my $exon_str = join(":", map{$_->start."-".$_->end} @{$tr_o->get_all_Exons});
#     $old_exon_str{$exon_str}++;
#   }  
#   
#   my %new_trids;
#   my %new_exon_str;
#   foreach my $tr_n (@{$gene_n->get_all_Transcripts}){
#     next if $tr_n->biotype =~ /^(artifact)$/;
#     my $tid = $tr_n->stable_id;
#     unless ($old_trids{$tid}){
#       print OUT "NEW_TR1\t$gsymbol\t$ensid\t$ottid\t".get_ens_transcript_id($tid)."\t$tid\n";
#     }
#     my $exon_str = join(":", map{$_->start."-".$_->end} @{$tr_n->get_all_Exons});
#     unless ($old_exon_str{$exon_str}){
#       print OUT "NEW_TR2\t$gsymbol\t$ensid\t$ottid\t".get_ens_transcript_id($tid)."\t$tid\n";
#     }
#   }



#a) new exon
#- exon in new annotation that shares no sequence with any exon in old annotation - 
#  excluding retained_intron transcripts from the old_annotation (as in previous analysis) ?

  my %novel_exon_ids;
  my %trs_w_novel_feat;
  my @novel_exons;
  my %tr_novel_exon_status;
  my %tr_biotypes;
  my @novel_introns;
  my %gene_novel_feat;

  foreach my $tr_n (@{$gene_n->get_all_Transcripts}){
    next if $tr_n->biotype =~ /^(retained_intron|artifact)$/;
    my $new_exon;
    EXON:foreach my $ex_n (@{$tr_n->get_all_Exons}){
      $new_exon = 1;
      my $chr_name = $ex_n->seq_region_name;
      my $slice_o = $oa_sa->fetch_by_region("chromosome", $chr_name, $ex_n->start, $ex_n->end);
      #Get exons in old annotation overlapping the exon in new annotation
      foreach my $ex_o (grep {$_->strand==$ex_n->strand} @{$slice_o->get_all_Exons}){
        my @trs = @{$oa_ta->fetch_all_by_exon_stable_id($ex_o->stable_id)};
        foreach my $tr (@trs){
          #Skip exon in new annotation if it overlaps an exon in a non-retained-intron transcript in the old annotation
          if ($tr->get_Gene->stable_id eq $ensid and !($tr->biotype =~ /^(retained_intron)$/)){
            $new_exon = 0;
            next EXON;
          }
        }
      }
      if ($new_exon){
        push(@novel_exons, $ex_n);
        $novel_exon_ids{$ex_n->stable_id} = 1;
        $trs_w_novel_feat{new_exons}{$tr_n->stable_id}{$ex_n->start."-".$ex_n->end}++;
        $gene_novel_feat{new_exons}{$ex_n->start."-".$ex_n->end}++;
        print OUT $tr_n->stable_id."  ".$ex_n->stable_id."\t".$ex_n->length."\t";
        if ($tr_n->exon_rank($ex_n) == 1 or $tr_n->exon_rank($ex_n) == scalar(@{$tr_n->get_all_Exons})){
          print OUT "terminal exon\n";
          $tr_novel_exon_status{$tr_n->stable_id}{"terminal"}++;

        }
        else{
          print OUT "internal exon\n";
          $tr_novel_exon_status{$tr_n->stable_id}{"internal"}++;
        }
        print IDS4GTF get_ids_for_gtf($ex_n, $tr_n, 'novel_exon')."\n";
      }
    }
  }
  #instances
  print OUT "Number of novel exons: ".scalar(keys %novel_exon_ids)."\n";
  #genomic coverage (non-redundant)
  my $new_exon_coverage = get_nr_exon_coverage(\@novel_exons);
  print OUT "Coverage: ".$new_exon_coverage."\n";



#b) new introns
  #Novel introns - simple count of introns in new and not in v20
  #New analysis not done last time - can you also count novel exon-skipping events - introns in new that are a) missing in v20 and b) connect two splice sites that are present in v20 c) where the splice sites in b) have an exon between them in the transcripts where they are found together

  my %seen_introns; 
  foreach my $tr_o (@{$gene_o->get_all_Transcripts}){
    next if $tr_o->biotype =~ /^(retained_intron)$/;
    foreach my $intron_o (@{$tr_o->get_all_Introns}){
      $seen_introns{$intron_o->start."-".$intron_o->end} = 1;
    }
  }
  print OUT "Old introns: ".scalar(keys %seen_introns)." - ".join(",", keys %seen_introns)."\n";

  my %new_introns;
  my %exon_skipping_events;
  foreach my $tr_n (@{$gene_n->get_all_Transcripts}){
    next if $tr_n->biotype =~ /^(artifact)$/;
    foreach my $intron_n (@{$tr_n->get_all_Introns}){
      next if $seen_introns{$intron_n->start."-".$intron_n->end};
      $new_introns{$intron_n->start."-".$intron_n->end} = 1;
      push(@novel_introns, $intron_n);
      print IDS4GTF get_ids_for_gtf($intron_n, $tr_n, 'novel_intron_exon_skip')."\n";
      $trs_w_novel_feat{new_introns}{$tr_n->stable_id}{$intron_n->start."-".$intron_n->end}++;
      $gene_novel_feat{new_introns}{$intron_n->start."-".$intron_n->end}++;      
      
      foreach my $tr_o (@{$gene_o->get_all_Transcripts}){
        #next if $tr_o->biotype =~ /^(retained_intron)$/;
        my ($start_ok, $end_ok, $seen);
        foreach my $intron_o (@{$tr_o->get_all_Introns}){
          if ($intron_o->start==$intron_n->start and !($intron_o->end==$intron_n->end)){
            $start_ok = 1;
          }
          elsif ($intron_o->end==$intron_n->end and !($intron_o->start==$intron_n->start)){
            $end_ok = 1;
          }
          elsif ($intron_o->start==$intron_n->start and $intron_o->end==$intron_n->end){
            $seen = 1;
          }
        }
        #Report if intron start and end are present in old annotation but not in the same intron
	#NOTE: this will report transcripts where one or multiple exons have been skipped
        if ($start_ok and $end_ok and !$seen){
          #print OUT "New exon skip: ".$intron_n->start."-".$intron_n->end." ".$tr_n->stable_id." ".$tr_n->biotype."\n";
          $exon_skipping_events{$intron_n->start."-".$intron_n->end} = 1;
          #push(@novel_introns, $intron_n);
          print IDS4GTF get_ids_for_gtf($intron_n, $tr_n, 'exon_skip')."\n";
          $trs_w_novel_feat{exon_skip}{$tr_n->stable_id}{$intron_n->start."-".$intron_n->end}++;
          $gene_novel_feat{exon_skip}{$intron_n->start."-".$intron_n->end}++;
        }
      }
    }
  }

#count:
#instances
  print OUT "New introns: ".scalar(keys %new_introns)." - ".join(",", sort keys %new_introns)."\n";
  print OUT "New introns (exon skips): ".scalar(keys %exon_skipping_events)." - ".join(",", sort keys %exon_skipping_events)."\n";




#c) shifted splice junction

#- exon that shares some (>1 base) sequence overlap with an exon in old annotation but has a new splice site - NB must be new splice site (i.e. connects to upstream and/or downstream exon) and not simply new exon boundary - NB2 please exclude retained_intron transcripts from the old_annotation from this analysis
#Count for shifted splice junctions (exons in new sharing overlap with v20 annotation but have different splice sites) - need to be sure to count splice sites (i.e. with connections to other exons) and exclude partial terminal exons - if an exon has two novel splice sites compared to overlapping v20 exon count 2

  my %seen_sjs_3; 
  my %seen_introns_3;
  foreach my $tr_o (@{$gene_o->get_all_Transcripts}){
    next if $tr_o->biotype =~ /^(retained_intron)$/;
    foreach my $intron_o (@{$tr_o->get_all_Introns}){
      $seen_sjs_3{$intron_o->start-1} = 1;
      $seen_sjs_3{$intron_o->end+1} = 1;
      $seen_introns_3{($intron_o->start-1)."-".($intron_o->end+1)} = 1;
    }
  }
  my %shifted_sjs_3;
  my @shifted_sj_feats_3;
  foreach my $tr_n (@{$gene_n->get_all_Transcripts}){
    next if $tr_n->biotype =~ /^(retained_intron|artifact)$/;

    foreach my $intron_n (@{$tr_n->get_all_Introns}){
      next if $seen_introns_3{($intron_n->start-1)."-".($intron_n->end+1)};
      my $prev_exon_n = $intron_n->prev_Exon;
      my $next_exon_n = $intron_n->next_Exon;
      
      foreach my $tr_o (@{$gene_o->get_all_Transcripts}){
        next if $tr_o->biotype =~ /^(retained_intron)$/;
        foreach my $intron_o (@{$tr_o->get_all_Introns}){
	  my $extra_feat;
	  #Intron must overlap (but not be identical) between old and new annotation
          if ($intron_o->start<=$intron_n->end and $intron_o->end>=$intron_n->start and !($intron_o->start==$intron_n->start and $intron_o->end==$intron_n->end)){
            my $prev_exon_o = $intron_o->prev_Exon;
            my $next_exon_o = $intron_o->next_Exon;
	
	    #Check intron 5' end
	    if ($prev_exon_o->start<=$prev_exon_n->end and $prev_exon_o->end>=$prev_exon_n->start){
              if ($tr_n->strand==1){
	        unless ($seen_sjs_3{$prev_exon_n->end}){
	          $shifted_sjs_3{$prev_exon_n->end}++;                 
	          $trs_w_novel_feat{shifted_sj_3}{$tr_n->stable_id}{$prev_exon_n->end}++;
                  $gene_novel_feat{shifted_sj_3}{$prev_exon_n->end}++;  
	          print IDS4GTF get_ids_for_gtf($intron_n, $tr_n, 'shifted_sj')."\n";
	          if ($prev_exon_n->end > $prev_exon_o->end){
		    #new: #####-------##### (+)
		    #old: ###---------##### (+)
                    $extra_feat = new Bio::EnsEMBL::Feature( -start => $prev_exon_o->end+1, -end => $prev_exon_n->end, -strand => 1, -slice => $prev_exon_n->slice, -analysis =>  $prev_exon_n->analysis ) or die "can't create new feature";
                  }
	        }
	      }
              elsif ($tr_n->strand==-1){
	        unless ($seen_sjs_3{$prev_exon_n->start}){
	          $shifted_sjs_3{$prev_exon_n->start}++;                 
	          $trs_w_novel_feat{shifted_sj_3}{$tr_n->stable_id}{$prev_exon_n->start}++;
                  $gene_novel_feat{shifted_sj_3}{$prev_exon_n->start}++; 
		  print IDS4GTF get_ids_for_gtf($intron_n, $tr_n, 'shifted_sj')."\n";
		  if ($prev_exon_n->start < $prev_exon_o->start){
		    #new: #####-------##### (-)
		    #old: #####---------### (-)		  		  
                    $extra_feat = new Bio::EnsEMBL::Feature( -start => $prev_exon_n->start, -end => $prev_exon_o->start-1, -strand => -1, -slice => $prev_exon_n->slice, -analysis =>  $prev_exon_n->analysis ) or die "can't create new feature";
                  }	      
	        }
	      }	      
	    }
	    #Check intron 3' end
            if ($next_exon_o->start<=$next_exon_n->end and $next_exon_o->end>=$next_exon_n->start){
              if ($tr_n->strand==1){
	        unless ($seen_sjs_3{$next_exon_n->start}){
	          $shifted_sjs_3{$next_exon_n->start}++;
	          $trs_w_novel_feat{shifted_sj_3}{$tr_n->stable_id}{$next_exon_n->start}++;
                  $gene_novel_feat{shifted_sj_3}{$next_exon_n->start}++;
		  print IDS4GTF get_ids_for_gtf($intron_n, $tr_n, 'shifted_sj')."\n";
	          if ($next_exon_n->start < $next_exon_o->start){
		    #new: #####-------##### (+)
		    #old: #####---------### (+)
                    $extra_feat = new Bio::EnsEMBL::Feature( -start => $next_exon_n->start, -end => $next_exon_n->start-1, -strand => 1, -slice => $next_exon_n->slice, -analysis =>  $next_exon_n->analysis ) or die "can't create new feature";
                  }		  		  
	        }
	      }
              elsif ($tr_n->strand==-1){
	        unless ($seen_sjs_3{$next_exon_n->end}){
	          $shifted_sjs_3{$next_exon_n->end}++;
	          $trs_w_novel_feat{shifted_sj_3}{$tr_n->stable_id}{$next_exon_n->end}++;
                  $gene_novel_feat{shifted_sj_3}{$next_exon_n->end}++;
		  print IDS4GTF get_ids_for_gtf($intron_n, $tr_n, 'shifted_sj')."\n";
	          if ($next_exon_n->end > $next_exon_o->end){
		    #new: #####-------##### (-)
		    #old: ###---------##### (-)
                    $extra_feat = new Bio::EnsEMBL::Feature( -start => $next_exon_o->end+1, -end => $next_exon_n->end, -strand => -1, -slice => $next_exon_n->slice, -analysis =>  $next_exon_n->analysis ) or die "can't create new feature";
                  }		  		  	      	       
	        }
              }
            }
	    push (@shifted_sj_feats_3, $extra_feat) if $extra_feat;
          }
        }
      }
    }
  }



###



#d) new CDS coverage - how many CDS genomic positions in new are not CDS in old (it can be from non-coding in old to CDS in new)
  my $cds_cov_n = get_cds_cov($gene_n);
  my $cds_cov_o = get_cds_cov($gene_o);
  my $new_cds_coverage = 0;
  foreach my $pos (sort keys %$cds_cov_n){
    $new_cds_coverage++ unless $cds_cov_o->{$pos};
  }
  
#Number of extra bases in new annotation, classified as exon, CDS, 5' UTR and 3' UTR
#Unlike the above one, the CDS section counts how many of the newly covered genomic positions correspond to a CDS - this is more restrictive than the above
  my %extra_exon_bases;
  my %extra_cds_bases;
  my %extra_5utr_bases;
  my %extra_3utr_bases;  
  my %cov_o;
  foreach my $tr_o (@{$gene_o->get_all_Transcripts}){
    #next if $tr_o->biotype =~ /^(retained_intron|artifact)$/;   
    foreach my $ex_o (@{$tr_o->get_all_Exons}){
      for (my $i=$ex_o->start; $i<=$ex_o->end; $i++){
        $cov_o{$i} = 1;
      }
    }
  }     
  #Exon 
  foreach my $tr_n (@{$gene_n->get_all_Transcripts}){
    next if $tr_n->biotype =~ /^(retained_intron|artifact)$/;
    foreach my $ex_n (@{$tr_n->get_all_Exons}){
      for (my $i=$ex_n->start; $i<=$ex_n->end; $i++){
        unless ($cov_o{$i}){
	  $extra_exon_bases{$i} = 1;	
	}
      }
    }
    #5' UTR
    foreach my $utr5 (@{$tr_n->get_all_five_prime_UTRs}){
      for (my $i=$utr5->start; $i<=$utr5->end; $i++){
        if ($extra_exon_bases{$i}){
	  $extra_5utr_bases{$i} = 1;
        }
      }
    }
    #3' UTR
    foreach my $utr3 (@{$tr_n->get_all_three_prime_UTRs}){
      for (my $i=$utr3->start; $i<=$utr3->end; $i++){
        if ($extra_exon_bases{$i}){
	  $extra_3utr_bases{$i} = 1;
        }
      }
    }              
  }
     
  #CDS
  foreach my $tr_n (@{$gene_n->get_all_Transcripts}){
    next if $tr_n->biotype =~ /^(retained_intron|artifact)$/;
    foreach my $exon (@{$tr_n->get_all_translateable_Exons}){
      for (my $i=$exon->start; $i<=$exon->end; $i++){
        if ($extra_exon_bases{$i}){
	  $extra_cds_bases{$i} = 1;
	  if ($extra_5utr_bases{$i}){
	    delete($extra_5utr_bases{$i});
	  }
	  if ($extra_3utr_bases{$i}){
	    delete($extra_3utr_bases{$i});
	  }	  
	}
      }
    }
  }




#count:
#instances 
  print OUT "Shifted splice junctions: ".scalar(keys %shifted_sjs_3)." - ".join(",", sort keys %shifted_sjs_3)."\n";    
  my $shifted_sj_extra_coverage = get_nr_exon_coverage(\@shifted_sj_feats_3);
  print OUT "Coverage: ".$shifted_sj_extra_coverage."\n";


  #Gene stats
  print OUT "GENE\t".$gsymbol."\t".$ensid."\t".
            scalar (keys %{$gene_novel_feat{new_exons}})."\t".
            scalar (keys %{$gene_novel_feat{new_introns}})."\t".
            scalar (keys %{$gene_novel_feat{exon_skip}})."\t".	    
	    scalar (keys %{$gene_novel_feat{shifted_sj_3}})."\t".	    
            $new_exon_coverage."\t".$shifted_sj_extra_coverage."\t".
            $new_cds_coverage."\t".
	    scalar (keys %extra_exon_bases)."\t".
	    scalar (keys %extra_cds_bases)."\t".
	    scalar (keys %extra_5utr_bases)."\t".
	    scalar (keys %extra_3utr_bases).
	    "\n";


  #Transcript stats
  foreach my $tr_n (@{$gene_n->get_all_Transcripts}){
    next if $tr_n->biotype =~ /^(artifact|retained_intron)$/;
    my $ens_tid = $tr_n->stable_id;
    print OUT "TR\t".$gsymbol."\t".$ensid."\t".$ens_tid."\t".$tr_n->biotype."\t".
              join(",", map {$_->code} grep {$_->code =~ /start_NF|end_NF/ and $_->value==1} @{$tr_n->get_all_Attributes()})."\t".
              join(",", keys %{$tr_novel_exon_status{$tr_n->stable_id}})."\t".
	      get_created_date($tr_n->stable_id)."\t".
              scalar(keys %{$trs_w_novel_feat{new_exons}{$tr_n->stable_id}})."\t".
              scalar(keys %{$trs_w_novel_feat{new_introns}{$tr_n->stable_id}})."\t".
              scalar(keys %{$trs_w_novel_feat{exon_skip}{$tr_n->stable_id}})."\t". 
	      scalar(keys %{$trs_w_novel_feat{shifted_sj_3}{$tr_n->stable_id}}).
	      "\n";
  }




}
close (OUT);
close (IDS4GTF);
close (LIST);


#########################

sub get_nr_exon_coverage {
  my $exons = shift; #set of non-unique exons
  my %cov;
  foreach my $exon (@$exons){
    for (my $i=$exon->start; $i<=$exon->end; $i++){
      $cov{$i} = 1;
    }
  }
  return scalar(keys %cov);
}


sub get_top_biotype {
  my $biotypes = shift;
  if ($biotypes->{"protein_coding"}){
    return "protein_coding";
  }
  elsif ($biotypes->{"nonsense_mediated_decay"}){
    return "nonsense_mediated_decay";
  }
  elsif ($biotypes->{"processed_transcript"}){
    return "processed_transcript";
  }
  elsif ($biotypes->{"retained_intron"}){
    return "retained_intron";
  }
  return 0;
}


###
###Print GTF for the whole gene list and grep exon lines for exon id + transcript id; exon ids of novel and flanking exons must be obtained previously
sub get_ids_for_gtf {
  my ($exon, $tr, $type) = @_;
  my @lines;
  if ($type eq 'novel_exon'){
    my $rank = $tr->exon_rank($exon);
    foreach my $ex (@{$tr->get_all_Exons}){
      if ($tr->exon_rank($ex) == $rank-1 or $tr->exon_rank($ex) == $rank or $tr->exon_rank($ex) == $rank+1){
        push(@lines, $type."\t".$tr->stable_id."\t".$ex->stable_id."\t");
      }
    }
  }
  elsif ($type eq 'novel_intron' or $type eq 'shifted_sj'){
    my $intron = $exon;
    push(@lines, $type."\t".$tr->stable_id."\t".$intron->prev_Exon->stable_id."\t");
    push(@lines, $type."\t".$tr->stable_id."\t".$intron->next_Exon->stable_id."\t");
  }
  return join("\n", @lines);
}


sub get_cds_cov {
  my $gene = shift;
  my %cds_cov;
  foreach my $tr (@{$gene ->get_all_Transcripts}){
    next if $tr->biotype =~ /^(retained_intron|artifact)$/;
    foreach my $exon (@{$tr->get_all_translateable_Exons}){
      for (my $i=$exon->start; $i<=$exon->end; $i++){
        $cds_cov{$i} = 1;
      }
    }
  }
  return \%cds_cov;
}


sub get_created_date {
  my $tid = shift;
  my $sql = "SELECT created_date FROM transcript t
             WHERE stable_id=? AND is_current";
  my $sth = $new_annot_db->dbc->prepare($sql);
  $sth->execute($tid);
  while (my ($datetime) = $sth->fetchrow_array){
    my ($date, $time) = split(/\s+/, $datetime);
    return $date;
  }
  return "NA";
}


