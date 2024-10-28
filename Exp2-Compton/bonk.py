# Import required modules
import struct

# Constants
HEADER_SIZE = 32
CHANNELS_PER_RECORD = 8

# Open the .CHN file (update the file path as necessary)
with open('cs_0.chn', 'rb') as f:
    # Read the first 32 bytes of header data
    header_data = f.read(HEADER_SIZE)
    
    # Unpack the header data (Fortran 'INTEGER*2' is 'int16' in Python and 'INTEGER*4' is 'int32')
    TYPE, MCA, SEG, SRTSEC, RLTIME, LVETME, SRTDTE, SRTTME, STRTCH, LNGTDT = struct.unpack('h h h 2s i i 8s 4s h h', header_data)
    
    # Check if the TYPE is -1 to verify if the file is valid
    if TYPE != -1:
        print("Invalid file type.")
    else:
        # Print header data
        print(f"TYPE = {TYPE}  MCA # = {MCA}  SEGMENT # = {SEG}")
        print(f"REALTIME = {RLTIME / 50:.2f} SECONDS, LIVETIME = {LVETME / 50:.2f} SECONDS")
        print(f"DATA COLLECTED AT {SRTTME.decode()}:{SRTSEC.decode()} ON {SRTDTE.decode()}")
        print(f"STARTING CHANNEL = {STRTCH}, NUMBER OF CHANNELS = {LNGTDT}")

        # Calculate the last record position
        LREC = 3 + (LNGTDT - 1) // CHANNELS_PER_RECORD

        # Seek to the first trailer record position
        f.seek(LREC * HEADER_SIZE)
        
        # Attempt to read 24 bytes for the trailer record
        trailer_record = f.read(24)  # Update this to the expected byte length (24)
        
        # Verify that we have read enough bytes before unpacking
        if len(trailer_record) == 24:
            TLRTYP, IS, ENG0, ENG1, X1, FW0, FW1 = struct.unpack('h h f f f f f', trailer_record)

            # Print energy calibration information
            print(f"ENERGY ZERO = {ENG0:.8f}")
            print(f"ENERGY SLOPE = {ENG1:.8f}")
            print(f"FWHM ZERO = {FW0:.8f}")
            print(f"FWHM SLOPE = {FW1:.8f}")
        else:
            print("Error: Could not read the full trailer record.")

        # Prompt for a channel number
        channel_num = int(input("Enter channel number: "))

        # Calculate block information for the specified channel
        channel = channel_num - 1
        BEGREC = channel // CHANNELS_PER_RECORD
        ENDREC = channel // CHANNELS_PER_RECORD

        # Loop through the records to read and display channel data
        for i in range(BEGREC + 2, ENDREC + 3):
            # Read a record of 8 channels
            f.seek((i - 1) * HEADER_SIZE)
            record_data = f.read(HEADER_SIZE)
            if len(record_data) == HEADER_SIZE:
                channels = struct.unpack('8h', record_data)  # Read 8 channels

                # Print the channels
                print(f"{1 + 8 * (i - 2):5}", *channels)
            else:
                print(f"Error: Could not read record at position {i}.")

# Close the file is automatically handled by 'with'
