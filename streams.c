#include "image_compression.h"

/*
int pushto_OUTSTREAM(OUTSTREAM *out, uint data, int bits) {
    // bits are the number of bits to encode from least to most significant
    if (bits < 1 || bits > 32) return 1;
    if (out->eof) return -1;
    data &= (bits == 32) ? ~0 : (1 << bits) - 1; // remove trash data, but if whole int is needed, we cant bitshift by 0
    int slip = bits + out->written_bits - 8;
    if (slip > 0) {
        out->buffer[out->written_bytes] |= data >> +slip;
        if (out->written_bytes >= out->buffer_bytes) { // overflow -> new buffer
            fwrite(out->buffer, 1, out->buffer_bytes, out->file);
            size_t els = 0; out->buffer_bytes++;
            while (els == 0 && out->buffer_bytes > 0) { // create the new biggest buffer possible
                out->buffer_bytes--;
                free(out->buffer); out->buffer = (char *)calloc(out->written_bytes, 1);
                els = fwrite(out->buffer, 1, out->buffer_bytes, out->file);
            } if (els == 0) { // if no buffer at all can be created, then we are at the end of the file
                out->eof = true;
            } else {
                out->written_bytes = -1; // if a buffer was created, reset the bytes counter (since im going to add 1)
            }
        }
        data &= (1 << slip) - 1; // remove alredy coded part
        out->written_bytes++; out->written_bits = 0;
        pushto_OUTSTREAM(out, data, slip);
    } else if (slip < 0) {
        out->buffer[out->written_bytes] |= data << -slip;
        out->written_bits += bits; 
    } else {
        out->buffer[out->written_bytes] |= data;
        out->written_bits = 0; out->written_bytes++;
    } 
    return 0;
}
*/
bool is_little_endian() {
    uint16_t a = 1;
    return *(char *)&a == 1;
}
void reset_bytes(OUTSTREAM *out) {
    if (out->written_bytes == (out->buffer_bytes-1)) { // there has been slip, and at we had alredy written up to the second to last byte
        fwrite(out->buffer, 1, out->buffer_bytes, out->file);
        out->written_bytes = -1;    
    }   out->written_bits = 0;
        out->written_bytes++;
}

int pushto_OUTSTREAM(OUTSTREAM *out, uint data, int bits) {
    // bits are the number of bits to encode from least to most significant
    if (bits < 1 || bits > 32) return 1;
    if (out->eof) return -1;
    data &= (bits == 32) ? ~0 : (1 << bits) - 1; // remove trash data, but if whole int is needed, we cant bitshift by 0
    int slip = bits + out->written_bits - 8;   
    if (slip > 0) {
        out->buffer[out->written_bytes] |= data >> +slip;
        reset_bytes(out);
        //data &= (1 << slip) - 1; // remove alredy coded part
        pushto_OUTSTREAM(out, data, slip);
    } else if (slip < 0) {
        out->buffer[out->written_bytes] |= data << -slip;
        out->written_bits += bits; 
    } else if (slip == 0) {
        out->buffer[out->written_bytes] |= data;
        reset_bytes(out);
    } 
    printf("Cache: \n");
    for (int i = 0; i < out->buffer_bytes; i++) {
        
        print_ubits(out->buffer[i]);
    }
    printf("\n");
    return 0;
}
OUTSTREAM *new_OUTSTREAM(const char* filename, int buffer_bytes) {
    if (buffer_bytes < 1) return NULL;
    OUTSTREAM *out = malloc(sizeof(OUTSTREAM));
    out->filename = filename;
    out->file = fopen(filename, "wb");
    out->written_bits = 0;
    out->written_bytes = 0;
    out->eof = false;
    out->buffer_bytes = buffer_bytes;
    out->buffer = calloc(buffer_bytes, 1); 
    //if (fread(out->buffer, 1, buffer_bytes, out->file)==0) {fclose(out->file); free(out); return NULL;}
    return out;
}

int delete_OUTSTREAM(OUTSTREAM *out) {
    if (out == NULL) return -1;
    if (out->written_bits != 0 | out->written_bytes!=0) fwrite(out->buffer, 1, out->written_bytes, out->file); // flush buffer if something has been written
    fclose(out->file);
    free(out->buffer);
    free(out);
    return 0;
}