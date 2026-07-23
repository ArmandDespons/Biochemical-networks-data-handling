from .Species import * 

__all__ = ["FERMENTATION_PRODUCTS"]

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