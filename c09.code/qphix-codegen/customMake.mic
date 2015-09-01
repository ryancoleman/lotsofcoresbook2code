mode=mic

# Generate AVX512 code for KNL/SKL
AVX512=0

# Define compute precision (1=float, 2=double)
PRECISION ?= 1

# Enable Lower Precision if set to 1
ENABLE_LOW_PRECISION ?= 0

# define SOA length 
# can be one of 4, 8 or 16 for MIC/SP
# and 2, 4 or 8 for MIC/DP
SOALEN ?= 8

# Enable serial spin compute in Dslash
SERIAL_SPIN=1

# Prefetching options
# FOR MIC SET THESE ALL TO 1
#
PREF_L1_SPINOR_IN = 1
PREF_L2_SPINOR_IN = 1
PREF_L1_SPINOR_OUT = 1
PREF_L2_SPINOR_OUT = 1
PREF_L1_GAUGE = 1
PREF_L2_GAUGE = 1
PREF_L1_CLOVER = 1
PREF_L2_CLOVER = 1

# Gather / Scatter options
USE_LDUNPK = 1            # Use loadunpack instead of gather
USE_PKST = 1              # Use packstore instead of scatter
USE_SHUFFLES = 0          # Use loads & Shuffles to transpose spinor when SOALEN>4
NO_GPREF_L1 = 1           # Generate bunch of normal prefetches instead of one gather prefetch for L1
NO_GPREF_L2 = 1           # Generate bunch of normal prefetches instead of one gather prefetch for L2

# Enable nontemporal streaming stores
ENABLE_STREAMING_STORES ?= 1
USE_PACKED_GAUGES ?= 1     # Use 2D xy packing for Gauges
USE_PACKED_CLOVER ?= 1     # Use 2D xy packing for Clover
