/** THIS IS AN AUTOMATICALLY GENERATED FILE.  DO NOT MODIFY
 * BY HAND!!
 *
 * Generated by lcm-gen
 *
 * From Source File: /Users/ryan/dev/symforce-zig-build/lcmtypes/symforce.lcm
 **/

#include <lcm/lcm_coretypes.h>

#ifndef __sym_key_t_hpp__
#define __sym_key_t_hpp__

#include <ostream>

namespace sym
{

class key_t
{
    public:
        // TODO(hayk): Issue for alignment?
        uint8_t letter;

        int64_t subscript;

        int64_t superscript;

    public:
        key_t() = default;

        /**
         * Member constructor
         */
        inline key_t(
            const uint8_t& letter_arg,
            const int64_t& subscript_arg,
            const int64_t& superscript_arg
        );

        /**
         * Encode a message into binary form.
         *
         * @param buf The output buffer.
         * @param offset Encoding starts at thie byte offset into @p buf.
         * @param maxlen Maximum number of bytes to write.  This should generally be
         *  equal to getEncodedSize().
         * @return The number of bytes encoded, or <0 on error.
         */
        inline __lcm_buffer_size encode(void *buf, __lcm_buffer_size offset, __lcm_buffer_size maxlen) const;

        /**
         * Check how many bytes are required to encode this message.
         */
        inline __lcm_buffer_size getEncodedSize() const;

        /**
         * Decode a message from binary form into this instance.
         *
         * @param buf The buffer containing the encoded message.
         * @param offset The byte offset into @p buf where the encoded message starts.
         * @param maxlen The maximum number of bytes to read while decoding.
         * @return The number of bytes decoded, or <0 if an error occured.
         */
        inline __lcm_buffer_size decode(const void *buf, __lcm_buffer_size offset, __lcm_buffer_size maxlen);

        /**
         * Retrieve the 64-bit fingerprint identifying the structure of the message.
         * Note that the fingerprint is the same for all instances of the same
         * message type, and is a fingerprint on the message type definition, not on
         * the message contents.
         */
        constexpr static int64_t getHash();

        using type_name_array_t = const char[6];

        inline static constexpr type_name_array_t* getTypeNameArrayPtr();

        /**
         * Returns "key_t"
         */
        inline static constexpr const char* getTypeName();

        /**
         * Returns "sym"
         */
        inline static constexpr const char * getPackageName();

        // LCM support functions. Users should not call these
        inline __lcm_buffer_size _encodeNoHash(void *buf, __lcm_buffer_size offset, __lcm_buffer_size maxlen) const;
        inline __lcm_buffer_size _getEncodedSizeNoHash() const;
        inline __lcm_buffer_size _decodeNoHash(const void *buf, __lcm_buffer_size offset, __lcm_buffer_size maxlen);
        constexpr static uint64_t _computeHash(const __lcm_hash_ptr *)
        {
            uint64_t hash = 0xd4cf464177ae74aaLL;
            return (hash<<1) + ((hash>>63)&1);
        }

        // Comparison operators.
        inline bool operator==(const key_t& other) const;
        inline bool operator!=(const key_t& other) const;

        // Ability to print to standard streams as well as the fmt library.
        friend std::ostream& operator<<(std::ostream& stream, const key_t& obj) {
#if defined(SKYMARSHAL_PRINTING_ENABLED)
            stream << "key_t(";
            stream << "letter=" << static_cast<int16_t>(obj.letter) << ", ";
            stream << "subscript=" << obj.subscript << ", ";
            stream << "superscript=" << obj.superscript;
            stream << ")";
#else
            stream << "<FORMATTING DISABLED>";
#endif
            return stream;
        }
};

key_t::key_t(
    const uint8_t& letter_arg,
    const int64_t& subscript_arg,
    const int64_t& superscript_arg
) : letter(letter_arg),
    subscript(subscript_arg),
    superscript(superscript_arg) {}

__lcm_buffer_size key_t::encode(void *buf, __lcm_buffer_size offset, __lcm_buffer_size maxlen) const
{
    __lcm_buffer_size pos = 0, tlen;
    int64_t hash = (int64_t)getHash();

    tlen = __int64_t_encode_array(buf, offset + pos, maxlen - pos, &hash, 1);
    if(tlen < 0) return tlen; else pos += tlen;

    tlen = this->_encodeNoHash(buf, offset + pos, maxlen - pos);
    if (tlen < 0) return tlen; else pos += tlen;

    return pos;
}

__lcm_buffer_size key_t::decode(const void *buf, __lcm_buffer_size offset, __lcm_buffer_size maxlen)
{
    __lcm_buffer_size pos = 0, thislen;

    int64_t msg_hash;
    thislen = __int64_t_decode_array(buf, offset + pos, maxlen - pos, &msg_hash, 1);
    if (thislen < 0) return thislen; else pos += thislen;
    if (msg_hash != getHash()) return -1;

    thislen = this->_decodeNoHash(buf, offset + pos, maxlen - pos);
    if (thislen < 0) return thislen; else pos += thislen;

    return pos;
}

__lcm_buffer_size key_t::getEncodedSize() const
{
    return 8 + _getEncodedSizeNoHash();
}

constexpr int64_t key_t::getHash()
{
    return static_cast<int64_t>(_computeHash(NULL));
}

constexpr key_t::type_name_array_t* key_t::getTypeNameArrayPtr() {
    return &"key_t";
}

constexpr const char* key_t::getTypeName()
{
    return *key_t::getTypeNameArrayPtr();
}

constexpr const char * key_t::getPackageName()
{
    return "sym";
}

__lcm_buffer_size key_t::_encodeNoHash(void *buf, __lcm_buffer_size offset, __lcm_buffer_size maxlen) const
{
    __lcm_buffer_size pos = 0, tlen;

    tlen = __byte_encode_array(buf, offset + pos, maxlen - pos, &this->letter, 1);
    if(tlen < 0) return tlen; else pos += tlen;

    tlen = __int64_t_encode_array(buf, offset + pos, maxlen - pos, &this->subscript, 1);
    if(tlen < 0) return tlen; else pos += tlen;

    tlen = __int64_t_encode_array(buf, offset + pos, maxlen - pos, &this->superscript, 1);
    if(tlen < 0) return tlen; else pos += tlen;

    return pos;
}

__lcm_buffer_size key_t::_decodeNoHash(const void *buf, __lcm_buffer_size offset, __lcm_buffer_size maxlen)
{
    __lcm_buffer_size pos = 0, tlen;

    tlen = __byte_decode_array(buf, offset + pos, maxlen - pos, &this->letter, 1);
    if(tlen < 0) return tlen; else pos += tlen;

    tlen = __int64_t_decode_array(buf, offset + pos, maxlen - pos, &this->subscript, 1);
    if(tlen < 0) return tlen; else pos += tlen;

    tlen = __int64_t_decode_array(buf, offset + pos, maxlen - pos, &this->superscript, 1);
    if(tlen < 0) return tlen; else pos += tlen;

    return pos;
}

__lcm_buffer_size key_t::_getEncodedSizeNoHash() const
{
    __lcm_buffer_size enc_size = 0;
    enc_size += __byte_encoded_array_size(NULL, 1);
    enc_size += __int64_t_encoded_array_size(NULL, 1);
    enc_size += __int64_t_encoded_array_size(NULL, 1);
    return enc_size;
}

bool key_t::operator==(const key_t& other) const {
  return (
          (letter==other.letter) && 
          (subscript==other.subscript) && 
          (superscript==other.superscript));
}

bool key_t::operator!=(const key_t& other) const {
  return !(*this==other);
}

}

#endif