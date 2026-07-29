from .Species import * 

__all__ = [
    "FERMENTATION_PRODUCTS",
    "FERMENTATION_CLASS"
]

FERMENTATION_PRODUCTS = Species.from_list([
    'acetate', 
    'formate', 
    'ethanol', 
    'lactate',
    'glycerol', 
    'pyruvate', 
    'butyrate', 
    'succinate', 
    ])

FERMENTATION_CLASS = {
    # =========================================================================
    # CLASS I: Non-oxidative Pyruvate Cleavage
    # (Junction yields 0 net/released electrons; oxidation stops at pyruvate)
    # =========================================================================

    # Yeasts (Ethanologens via PDC)
    "Saccharomyces": "I",
    "Schizosaccharomyces": "I",
    "Kluyveromyces": "I",
    "Pichia": "I",
    
    # Enteric & Mixed-Acid Fermenters (via PFL/LDH)
    "Escherichia": "I",
    "Salmonella": "I",
    "Klebsiella": "I",
    
    # Lactic Acid Bacteria (via LDH or PFL/Phosphoketolase)
    "Lactococcus": "I",
    "Streptococcus": "I",
    "Lactobacillus": "I",
    "Leuconostoc": "I",
    "Pediococcus": "I",
    "Enterococcus": "I",
    "Bifidobacterium": "I",
    

    # =========================================================================
    # CLASS II: oxidative Pyruvate Cleavage
    # (Junction actively liberates electrons onto Ferredoxin/carriers via PFOR)
    # =========================================================================

    # Classic Clostridia & Acetogens
    "Clostridium": "II",
    "Acetobacterium": "II",
    
    # Core Gut & Rumen Anaerobes
    "Bacteroides": "II",          # Human gut dominant; utilizes PFOR to generate Acetyl-CoA
    "Ruminococcus": "II",         # Major rumen cellulolytic genus; relies heavily on PFOR 
    "Megasphaera": "II",          # Rumen/gut lactate-utilizer (e.g., M. elsdenii); PFOR driven
    "Veillonella": "II",          # Lactate-fermenting structural anaerobe; utilizes PFOR
    
    # Thermophilic Anaerobes
    "Thermoanaerobacter": "II",   # Highly active PFOR
    "Thermoanaerobacterium": "II" # Closely related pentose/hexose thermophilic fermenter
}