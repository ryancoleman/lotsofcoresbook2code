"""These functions can be helpful to have attached to
an atoms object.
"""
import types


def get_raw_score(self):
    return self.info['key_value_pairs']['raw_score']
    
    
def set_raw_score(self, score):
    self.info['key_value_pairs']['raw_score'] = score
    
    
def enable_raw_score_methods(a):
    if not 'key_value_pairs' in a.info:
        a.info['key_value_pairs'] = {}
    a.set_raw_score = types.MethodType(set_raw_score, a)
    a.get_raw_score = types.MethodType(get_raw_score, a)
