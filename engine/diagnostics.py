"""Cheap diagnostics used before expensive coupled campaigns."""
def friedmann_constraint(H,rho):
    return 3.0*H*H-rho

def exchange_closure(recipient,donor):
    return recipient+donor

def relative_error(a,b,floor=1e-30):
    return abs(a-b)/max(abs(a),abs(b),floor)
