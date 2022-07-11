#include "image_compression.h"

int pushto_OUTSTREAM(OUTSTREAM *out, uint data, int bits) {
    // bits are the number of bits to encode from least to most significant
    if (bits < 1 || bits > 32) return 1;
    if (out->eof) return -1;
    data &= (bits == 32) ? ~0 : (1 << bits) - 1; // remove trash data, but if whole int is needed, we cant bitshift by 0
    int slip = bits + out->d_bits - 8;
    if (slip > 0) {
        out->buffer[out->d_bytes] |= data >> +slip;
        if (out->d_bytes >= out->buffer_bytes) { // overflow -> new buffer
            fwrite(out->buffer, 1, out->buffer_bytes, out->file);
            size_t els = 0; out->buffer_bytes++;
            while (els == 0 && out->buffer_bytes > 0) { // create the new biggest buffer possible
                out->buffer_bytes--;
                free(out->buffer); out->buffer = (char *)calloc(out->d_bytes, 1);
                els = fwrite(out->buffer, 1, out->buffer_bytes, out->file);
            } if (els == 0) { // if no buffer at all can be created, then we are at the end of the file
                out->eof = true;
            } else {
                out->d_bytes = -1; // if a buffer was created, reset the bytes counter (since im going to add 1)
            }
        }
        data &= (1 << slip) - 1; // remove alredy coded part
        out->d_bytes++; out->d_bits = 0;
        pushto_OUTSTREAM(out, data, slip);
    } else if (slip < 0) {
        out->buffer[out->d_bytes] |= data << -slip;
        out->d_bits += bits; 
    } else {
        out->buffer[out->d_bytes] |= data;
        out->d_bits = 0; out->d_bytes++;
    } 
    return 0;
}

OUTSTREAM *new_OUTSTREAM(const char* filename, int buffer_bytes) {
    if (buffer_bytes < 1) return NULL;
    OUTSTREAM *out = (OUTSTREAM *)malloc(sizeof(OUTSTREAM));
    out->filename = filename;
    out->file = fopen(filename, "wb");
    out->d_bits = 0;
    out->d_bytes = 0;
    out->eof = false;
    out->buffer_bytes = buffer_bytes;
    out->buffer = (char *)calloc(buffer_bytes, 1); 
    //if (fread(out->buffer, 1, buffer_bytes, out->file)==0) {fclose(out->file); free(out); return NULL;}
    return out;
}

int delete_OUTSTREAM(OUTSTREAM *out) {
    if (out == NULL) return -1;
    if (out->d_bits != 0 | out->d_bytes!=0) fwrite(out->buffer, 1, out->buffer_bytes, out->file); // flush buffer if something has been written
    fclose(out->file);
    free(out->buffer);
    free(out);
    return 0;
}