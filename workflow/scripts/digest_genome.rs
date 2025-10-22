use clap::Parser;
use std::io::{Write,BufWriter};
use std::path::PathBuf;
use bio::io::fasta;
use bio::io::fasta::{Reader,Record};
use regex::Regex;
use std::error::Error;
use std::fs::File;

#[derive(Parser, Debug)]
#[clap(author = "GilbertHan", version, about = "Finds restriction fragments in a FASTA file")]
struct Cli {
    /// Input FASTA file
    #[clap(value_parser)]
    fastafile: PathBuf,

    /// One or more restriction sites (e.g., A^AGCTT)
    #[clap(short, long, required = true, num_args(1..))]
    restriction_sites: Vec<String>,

    /// Output BED file [default: <fastafile_prefix>_fragments.bed]
    #[clap(short, long)]
    out: Option<PathBuf>,
}

fn parse_restriction_sites(raw_site : &[String]) -> (Vec<String>, Vec<usize>){
    let mut all_site = Vec::new();
    let mut all_offsite = Vec::new();

    for site_str in raw_site.iter().flat_map(|s| s.split(',')){
        let offsite = site_str.find('^').unwrap_or_else(||{
            panic!("The input string {} must have '^'", site_str);
        });
        let sequence: String = site_str.chars().filter(|&c| c!= '^').collect();
        let expanded_sequences = expand_n(&sequence);
        let num_expanded = expanded_sequences.len();
        all_site.extend(expanded_sequences);
        all_offsite.extend(std::iter::repeat(offsite).take(num_expanded));
    }
    (all_site, all_offsite)
}

fn expand_n(sequence: &str) -> Vec<String> {
    if let Some(n_pos) = sequence.find('N') {
        let mut results = Vec::new();
        let (prefix, suffix) = sequence.split_at(n_pos);
        let suffix_no_n = &suffix[1..]; // Skip the 'N'

        for nuc in ['A', 'C', 'G', 'T'] {
            let new_seq = format!("{}{}{}", prefix, nuc, suffix_no_n);
            // Recurse to handle multiple 'N's
            results.extend(expand_n(&new_seq)); 
        }
        results
    } else {
        // Base case: no 'N' found
        vec![sequence.to_string()]
    }
}


fn main() -> Result<(), Box<dyn Error>> {
    let cli = Cli::parse();
    let (site_list, offset_list) = parse_restriction_sites(&cli.restriction_sites);
    // Regex pattern
    let patterns:Vec<Regex> = site_list.iter()
        .map(|seq|{
            let pattern_str = format!("(?i){}",seq);
            Regex::new(&pattern_str).expect("Fail to compile")
        })
        .collect();

    let file_path = cli.fastafile;
    let out_path = cli.out.unwrap_or_else(|| {
        let stem = file_path.file_stem().unwrap().to_str().unwrap_or("fasta");
        let mut p = file_path.clone();
        p.set_file_name(format!("{}_fragments.bed", stem));
        p
    });

    let fastaReader = Reader::from_file(&file_path)?;
    let file = File::create(out_path)?;
    let mut writer = BufWriter::new(file);
    let mut record = Record::new();
    let rec = fastaReader.records();
    for result in rec {
        let record = result?;
        let chrName = record.id();
        let seq = record.seq();
        let seq_len = seq.len();
        let seq_str = String::from_utf8_lossy(seq);
        let mut cut_site: Vec<usize> = Vec::new();
        for (i,re) in patterns.iter().enumerate(){
            let offset = offset_list[i];
            let mut start_index = 0;
            while let Some(mat) = re.find_at(&seq_str, start_index) {
                cut_site.push(mat.start() + offset);
                start_index = mat.start() + 1;
            }
        }
        cut_site.sort_unstable();
        let mut frag_id = 0;
        let mut last_pos = 0;

        for site in cut_site.iter(){
            let start = last_pos;
            let end = *site;
            if end > start {
                frag_id += 1;
                let frag_name = format!("HiC_{}_{}", chrName, frag_id);
                writeln!(writer,"{}\t{}\t{}\t{}\t0\t+", chrName, start, end, frag_name)?;
            }
            last_pos = *site;
        }

        if seq_len > last_pos {
            frag_id += 1;
            let frag_name = format!("HiC_{}_{}", chrName, frag_id);
            writeln!(writer,"{}\t{}\t{}\t{}\t0\t+", chrName, last_pos, seq_len, frag_name)?;
        }

    }

    Ok(())

}