#!/software/bin/perl -w

#This script compares the variants falling on annotation from two different releases (restricted to a list of genes)

# use the Ensembl dataset “All phenotype-associated - short variants (SNPs and indels)”

#1 - compare variants that intersect with pre-existing annotation and novel annotation
#2 - do a second pass to look at features that lie within 8 bases of flanking intronic sequence of old and new annotation providing they are not picked up in 1

#Report:
#1 - overall numbers -eg 1500 variants overlap features in preexisting annotation and 50 overlap new annotation
#2 - numbers broken down by source database for old/new annotation- e.g. clinvar, LSDB, dbSNP - this can be redundant as many variants are in multiple databases
#3 - a list of variants that overlap only novel features  - SNP IDs (multiple if they are in multiple datasets) and chr, strand, co-ordinates

use strict;
use Getopt::Long;
use Bio::EnsEMBL::DBSQL::DBAdaptor;;
use Bio::EnsEMBL::Registry;

$|=1;

my $infile;
my $outfile;
my $species;
my $ext = 8;

&GetOptions(
            'in=s'      => \$infile,
            'out=s'     => \$outfile,
            'species=s' => \$species,
            'ext=s'     => \$ext
           );

my $registry = 'Bio::EnsEMBL::Registry';
$registry->load_registry_from_db(
  -host => 'ensembldb.ensembl.org',
  -user => 'anonymous',
);

my $vfa = $registry->get_adaptor('human', 'variation', 'variationfeature');
my $sa = $registry->get_adaptor('human', 'core', 'slice');


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
my $oa_ga = $old_annot_db->get_GeneAdaptor;

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
my $na_ga = $new_annot_db->get_GeneAdaptor;


open (OUT, ">$outfile") or die "Can't open $outfile\n";

print OUT "VF\tGene symbol\tVariation name\tVariation source\tChromosome\tStart\tEnd\tStrand\n";
print OUT "EVF\tGene symbol\tVariation name\tVariation source\tChromosome\tStart\tEnd\tStrand\n";
my %variants_o;
my %uniq_var_o;
my %variants_o_ext;
my %uniq_var_o_ext;
my %variants_n;
my %uniq_var_n;
my %variants_n_ext;
my %uniq_var_n_ext;
my %variants_nf;
my %uniq_var_nf;
my %variants_nf_ext;
my %uniq_var_nf_ext;
my %printed_vars;

open (LIST, $infile) or die "Can't open $infile\n";
while (<LIST>){
  chomp;
  my ($ensid, $gsymbol) = split(/\t/);
  my $gene_o = $oa_ga->fetch_by_stable_id($ensid) or die "Can't find old gene $ensid";
  my $gene_n = $na_ga->fetch_by_stable_id($ensid) or die "Can't find new gene $ensid";


  print OUT "\n\n#############\n\nChecking gene $ensid / $gsymbol\n";

  foreach my $exon_o (@{$gene_o->get_all_Exons}){
    my $exon_slice = $sa->fetch_by_region("chromosome", $exon_o->seq_region_name, $exon_o->seq_region_start, $exon_o->seq_region_end);
    foreach my $vf (@{$vfa->fetch_all_by_Slice($exon_slice)}){
      $variants_o{$vf->source_name}{$vf->variation_name} = $vf->seq_region_name.":".$vf->seq_region_start."-".$vf->seq_region_end;
      #$uniq_var_o{$vf->seq_region_name.":".$vf->seq_region_start."-".$vf->seq_region_end}++;
    }
    #8bp exon-flanking sequences
    my $exon_f1_slice = $sa->fetch_by_region("chromosome", $exon_o->seq_region_name, $exon_o->seq_region_start - $ext, $exon_o->seq_region_start - 1);
    my $exon_f2_slice = $sa->fetch_by_region("chromosome", $exon_o->seq_region_name, $exon_o->seq_region_end + 1, $exon_o->seq_region_end + $ext);   
    foreach my $vf (@{$vfa->fetch_all_by_Slice($exon_f1_slice)}, @{$vfa->fetch_all_by_Slice($exon_f2_slice)}){
      unless ($variants_o{$vf->source_name}{$vf->variation_name}){
        $variants_o_ext{$vf->source_name}{$vf->variation_name} = $vf->seq_region_name.":".$vf->seq_region_start."-".$vf->seq_region_end;
        #$uniq_var_o_ext{$vf->seq_region_name.":".$vf->seq_region_start."-".$vf->seq_region_end}++;
      }
    }
  }
  
  my %seen_exons;
  foreach my $tr_n (@{$gene_n->get_all_Transcripts}){
    next if $tr_n->biotype =~ /^(artifact)$/;
    foreach my $exon_n (@{$tr_n->get_all_Exons}){
      next if $seen_exons{$exon_n->stable_id};
      $seen_exons{$exon_n->stable_id} = 1;
      my $srname = $exon_n->seq_region_name;
      my $exon_slice = $sa->fetch_by_region("chromosome", $srname, $exon_n->seq_region_start, $exon_n->seq_region_end);
      foreach my $vf (@{$vfa->fetch_all_by_Slice($exon_slice)}){
        $variants_n{$vf->source_name}{$vf->variation_name} = $vf->seq_region_name.":".$vf->seq_region_start."-".$vf->seq_region_end;
        #$uniq_var_n{$vf->seq_region_name.":".$vf->seq_region_start."-".$vf->seq_region_end}++;

        unless ($variants_o{$vf->source_name}{$vf->variation_name}){
          $variants_nf{$vf->source_name}{$vf->variation_name} = $vf->seq_region_name.":".$vf->seq_region_start."-".$vf->seq_region_end;
          #$uniq_var_nf{$vf->seq_region_name.":".$vf->seq_region_start."-".$vf->seq_region_end}++;
          unless ($printed_vars{$gsymbol.":".$vf->source_name.":".$vf->variation_name}){
            print OUT "VF\t$gsymbol\t".$vf->variation_name."\t".$vf->source_name."\t".$vf->seq_region_name."\t".$vf->seq_region_start."\t".$vf->seq_region_end."\t".($exon_n->seq_region_strand ? "+":"-")."\n";
            print OUT "S\t".join(", ", @{$vf->variation->get_all_synonyms('ClinVar')})."\n";
            $printed_vars{$gsymbol.":".$vf->source_name.":".$vf->variation_name} = 1;
          }
        }
      }
      
      #8bp exon-flanking sequences
      my $exon_f1_slice = $sa->fetch_by_region("chromosome", $srname, $exon_n->seq_region_start - $ext, $exon_n->seq_region_start  - 1);   
      my $exon_f2_slice = $sa->fetch_by_region("chromosome", $srname, $exon_n->seq_region_end + 1, $exon_n->seq_region_end  + $ext);               
      foreach my $vf (@{$vfa->fetch_all_by_Slice($exon_f1_slice)}, @{$vfa->fetch_all_by_Slice($exon_f2_slice)}){
        unless ($variants_n{$vf->source_name}{$vf->variation_name}){
          $variants_n_ext{$vf->source_name}{$vf->variation_name} = $vf->seq_region_name.":".$vf->seq_region_start."-".$vf->seq_region_end;
          #$uniq_var_n_ext{$vf->seq_region_name.":".$vf->seq_region_start."-".$vf->seq_region_end}++;

          unless ($variants_o_ext{$vf->source_name}{$vf->variation_name} or $variants_nf{$vf->source_name}{$vf->variation_name}){
            $variants_nf_ext{$vf->source_name}{$vf->variation_name} = $vf->seq_region_name.":".$vf->seq_region_start."-".$vf->seq_region_end;
            #$uniq_var_nf{$vf->seq_region_name.":".$vf->seq_region_start."-".$vf->seq_region_end}++;
            unless ($printed_vars{$gsymbol.":".$vf->source_name.":".$vf->variation_name}){
              print OUT "EVF\t$gsymbol\t".$vf->variation_name."\t".$vf->source_name."\t".$vf->seq_region_name."\t".$vf->seq_region_start."\t".$vf->seq_region_end."\t".($exon_n->seq_region_strand ? "+":"-")."\n";
              $printed_vars{$gsymbol.":".$vf->source_name.":".$vf->variation_name} = 1;
            }
          }
          
        }
      }    
    }
  }
  
}
close (LIST);

print OUT "\n\n\nEXON SEQUENCES\n\n";
print OUT "Old annotation\n";
my $sum_o;
foreach my $source (sort keys %variants_o){
  print OUT $source."\t".scalar(keys %{$variants_o{$source}})."\n";
  $sum_o += scalar(keys %{$variants_o{$source}});
  foreach my $name (keys %{$variants_o{$source}}){
    $uniq_var_o{$variants_o{$source}{$name}} = 1;
  }
}
print OUT "total\t$sum_o\n";
print OUT "total unique\t".scalar(keys %uniq_var_o)."\n\n";

print OUT "New annotation\n";
my $sum_n;
foreach my $source (sort keys %variants_n){
  print OUT $source."\t".scalar(keys %{$variants_n{$source}})."\n";
  $sum_n += scalar(keys %{$variants_n{$source}});
  foreach my $name (keys %{$variants_n{$source}}){
    $uniq_var_n{$variants_n{$source}{$name}} = 1;
  }
}
print OUT "total\t$sum_n\n";
print OUT "total unique\t".scalar(keys %uniq_var_n)."\n\n";

print OUT "Novel features\n";
my $sum_nf;
foreach my $source (sort keys %variants_nf){
  print OUT $source."\t".scalar(keys %{$variants_nf{$source}})."\n";
  $sum_nf += scalar(keys %{$variants_nf{$source}});
  foreach my $name (keys %{$variants_nf{$source}}){
    $uniq_var_nf{$variants_nf{$source}{$name}} = 1;
  }  
}
print OUT "total\t$sum_nf\n";
print OUT "total unique\t".scalar(keys %uniq_var_nf)."\n\n";


print OUT "\n\n\nEXON FLANKING SEQUENCES ($ext NT)\n\n";
print OUT "Old annotation\n";
my $sum_o_ext;
foreach my $source (sort keys %variants_o_ext){
  print OUT $source."\t".scalar(keys %{$variants_o_ext{$source}})."\n";
  $sum_o_ext += scalar(keys %{$variants_o_ext{$source}});
  foreach my $name (keys %{$variants_o_ext{$source}}){
    $uniq_var_o_ext{$variants_o_ext{$source}{$name}} = 1;
  }
}
print OUT "total\t$sum_o_ext\n";
print OUT "total unique\t".scalar(keys %uniq_var_o_ext)."\n\n";

print OUT "New annotation\n";
my $sum_n_ext;
foreach my $source (sort keys %variants_n_ext){
  print OUT $source."\t".scalar(keys %{$variants_n_ext{$source}})."\n";
  $sum_n_ext += scalar(keys %{$variants_n_ext{$source}});
  foreach my $name (keys %{$variants_n_ext{$source}}){
    $uniq_var_n_ext{$variants_n_ext{$source}{$name}} = 1;
  }
}
print OUT "total\t$sum_n_ext\n";
print OUT "total unique\t".scalar(keys %uniq_var_n_ext)."\n\n";

print OUT "Novel features\n";
my $sum_nf_ext;
foreach my $source (sort keys %variants_nf_ext){
  print OUT $source."\t".scalar(keys %{$variants_nf_ext{$source}})."\n";
  $sum_nf_ext += scalar(keys %{$variants_nf_ext{$source}});
  foreach my $name (keys %{$variants_nf_ext{$source}}){
    $uniq_var_nf_ext{$variants_nf_ext{$source}{$name}} = 1;
  }
}
print OUT "total\t$sum_nf_ext\n";
print OUT "total unique\t".scalar(keys %uniq_var_nf_ext)."\n\n";



close (OUT);


#########################



