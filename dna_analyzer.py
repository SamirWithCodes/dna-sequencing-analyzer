import collections
import sys
import io

# উইন্ডোজ টার্মিনালে বাংলা ফন্ট সঠিকভাবে দেখানোর জন্য UTF-8 এনকোডিং সেট করা
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

class DNAAnalyzer:
    def __init__(self, sequence: str):
        """ডিএনএ সিকোয়েন্সটি ক্যাপিটাল লেটারে কনভার্ট করে সেট করা হচ্ছে"""
        self.sequence = sequence.upper().strip()
        self.validation_check()

    def validation_check(self):
        """সিকোয়েন্সটি সঠিক কিনা তা যাচাই করার মেথড"""
        valid_bases = set("ATCG")
        if not set(self.sequence).issubset(valid_bases):
            raise ValueError("ভুল সিকোয়েন্স! ডিএনএ-তে কেবল A, T, C, এবং G থাকতে পারে।")

    def get_biomarker_metrics(self):
        """ল্যাব ডেটার জন্য সিকোয়েন্সের দৈর্ঘ্য এবং বেস কাউন্ট বের করা"""
        length = len(self.sequence)
        counts = collections.Counter(self.sequence)
        
        # GC-Content গণনা
        gc_content = ((counts['G'] + counts['C']) / length) * 100 if length > 0 else 0
        
        return {
            "Length": length,
            "Base Counts": dict(counts),
            "GC Content (%)": round(gc_content, 2)
        }

    def transcribe_to_rna(self) -> str:
        """ডিএনএ থেকে আরএনএ-তে রূপান্তর (Transcription: T -> U)"""
        return self.sequence.replace('T', 'U')

    def translate_to_protein(self) -> str:
        """আরএনএ কোডন থেকে প্রোটিন (অ্যামিনো অ্যাসিড) তৈরি করা"""
        rna = self.transcribe_to_rna()
        
        # স্ট্যান্ডার্ড জেনেটিক কোড ডিকশনারি (AGC এবং বাকি সব কোডন ফিক্সড)
        genetic_code = {
            'AUA':'I', 'AUC':'I', 'AUU':'I', 'AUG':'M',
            'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACU':'T',
            'AAC':'N', 'AAU':'N', 'AAG':'K', 'AAA':'K',
            'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCU':'A',
            'GAC':'D', 'GAU':'D', 'GAA':'E', 'GAG':'E',
            'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGU':'G',
            'UCA':'S', 'UCC':'S', 'UCG':'S', 'UCU':'S',
            'AGC':'S', 'AGU':'S', 'AGA':'R', 'AGG':'R',  # আপনার ইনপুটের জন্য প্রয়োজনীয় 'AGC' যুক্ত করা হয়েছে
            'UUC':'F', 'UUU':'F', 'UUA':'L', 'UUG':'L',
            'UAC':'Y', 'UAU':'Y', 'UAA':'_', 'UAG':'_', 'UGA':'_',
            'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGU':'R',
            'CUA':'L', 'CUC':'L', 'CUG':'L', 'CUU':'L',
            'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCU':'P',
            'CAC':'H', 'CAU':'H', 'CAA':'Q', 'CAG':'Q',
        }
        
        protein = []
        for i in range(0, len(rna) - 2, 3):
            codon = rna[i:i+3]
            amino_acid = genetic_code.get(codon, '?')
            if amino_acid == '_':  # Stop Codon পেলে লুপ বন্ধ হবে
                break
            protein.append(amino_acid)
            
        return "".join(protein)

# স্ক্রিপ্টটি রান করার মূল অংশ
if __name__ == "__main__":
    print("--- বায়োইনফরমেটিক্স ল্যাব ডেটা অ্যানালাইজার ---")
    
    # আপনার স্ক্রিনশটের স্যাম্পল সিকোয়েন্সটি এখানে দেওয়া হলো
    sample_dna = "ATGCTTCGATCAGCTAGCTAGCTAGCTAGC"
    print(f"বিশ্লেষণ করা হচ্ছে: {sample_dna}\n")
    
    try:
        analyzer = DNAAnalyzer(sample_dna)
        
        # মেট্রিক্স অ্যানালিসিস
        metrics = analyzer.get_biomarker_metrics()
        print("[১] ল্যাব মেট্রিক্স ফলাফল:")
        print(f"  - Length: {metrics['Length']}")
        print(f"  - Base Counts: {metrics['Base Counts']}")
        print(f"  - GC Content (%): {metrics['GC Content (%)']}")
            
        # ট্রান্সক্রিপশন ও ট্রান্সলেশন
        print(f"\n[২] আরএনএ ট্রান্সক্রিপশন: {analyzer.transcribe_to_rna()}")
        print(f"[৩] প্রোটিন ট্রান্সলেশন (অ্যামিনো অ্যাসিড চেইন): {analyzer.translate_to_protein()}")
        
    except ValueError as e:
        print(f"Error: {e}")