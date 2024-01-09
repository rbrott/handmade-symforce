/** THIS IS AN AUTOMATICALLY GENERATED FILE.  DO NOT MODIFY
 * BY HAND!!
 *
 * Generated by lcm-gen
 *
 * From Source File: /Users/ryan/dev/symforce-zig-build/lcmtypes/symforce.lcm
 **/

#include <lcm/lcm_coretypes.h>

#ifndef __sym_imu_biases_t_hpp__
#define __sym_imu_biases_t_hpp__

#include <ostream>
#include "lcmtypes/eigen_lcm/Vector3d.hpp"

namespace sym
{

class imu_biases_t
{
    public:
        ::eigen_lcm::Vector3d accel_bias;

        ::eigen_lcm::Vector3d gyro_bias;

    public:
        imu_biases_t() = default;

        /**
         * Member constructor
         */
        inline imu_biases_t(
            const ::eigen_lcm::Vector3d& accel_bias_arg,
            const ::eigen_lcm::Vector3d& gyro_bias_arg
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
        inline static int64_t getHash();

        using type_name_array_t = const char[13];

        inline static constexpr type_name_array_t* getTypeNameArrayPtr();

        /**
         * Returns "imu_biases_t"
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
#if !defined(SKYDIO_DISABLE_LCM_NO_INLINE)
        __attribute__((noinline))
#endif
        static uint64_t _computeHash(const __lcm_hash_ptr *p)
        {
            const __lcm_hash_ptr *fp;
            for(fp = p; fp != NULL; fp = fp->parent)
                if(fp->v == imu_biases_t::getHash)
                    return 0;
            const __lcm_hash_ptr cp = { p, imu_biases_t::getHash };

            uint64_t hash = 0x176ba7d50765eac2LL +
                ::eigen_lcm::Vector3d::_computeHash(&cp) +
         ::eigen_lcm::Vector3d::_computeHash(&cp);

            return (hash<<1) + ((hash>>63)&1);
        }

        // Comparison operators.
        inline bool operator==(const imu_biases_t& other) const;
        inline bool operator!=(const imu_biases_t& other) const;

        // Ability to print to standard streams as well as the fmt library.
        friend std::ostream& operator<<(std::ostream& stream, const imu_biases_t& obj) {
#if defined(SKYMARSHAL_PRINTING_ENABLED)
            stream << "imu_biases_t(";
            stream << "accel_bias=" << obj.accel_bias << ", ";
            stream << "gyro_bias=" << obj.gyro_bias;
            stream << ")";
#else
            stream << "<FORMATTING DISABLED>";
#endif
            return stream;
        }
};

imu_biases_t::imu_biases_t(
    const ::eigen_lcm::Vector3d& accel_bias_arg,
    const ::eigen_lcm::Vector3d& gyro_bias_arg
) : accel_bias(accel_bias_arg),
    gyro_bias(gyro_bias_arg) {}

__lcm_buffer_size imu_biases_t::encode(void *buf, __lcm_buffer_size offset, __lcm_buffer_size maxlen) const
{
    __lcm_buffer_size pos = 0, tlen;
    int64_t hash = (int64_t)getHash();

    tlen = __int64_t_encode_array(buf, offset + pos, maxlen - pos, &hash, 1);
    if(tlen < 0) return tlen; else pos += tlen;

    tlen = this->_encodeNoHash(buf, offset + pos, maxlen - pos);
    if (tlen < 0) return tlen; else pos += tlen;

    return pos;
}

__lcm_buffer_size imu_biases_t::decode(const void *buf, __lcm_buffer_size offset, __lcm_buffer_size maxlen)
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

__lcm_buffer_size imu_biases_t::getEncodedSize() const
{
    return 8 + _getEncodedSizeNoHash();
}

int64_t imu_biases_t::getHash()
{
    static int64_t hash = _computeHash(NULL);
    return hash;
}

constexpr imu_biases_t::type_name_array_t* imu_biases_t::getTypeNameArrayPtr() {
    return &"imu_biases_t";
}

constexpr const char* imu_biases_t::getTypeName()
{
    return *imu_biases_t::getTypeNameArrayPtr();
}

constexpr const char * imu_biases_t::getPackageName()
{
    return "sym";
}

__lcm_buffer_size imu_biases_t::_encodeNoHash(void *buf, __lcm_buffer_size offset, __lcm_buffer_size maxlen) const
{
    __lcm_buffer_size pos = 0, tlen;

    tlen = this->accel_bias._encodeNoHash(buf, offset + pos, maxlen - pos);
    if(tlen < 0) return tlen; else pos += tlen;

    tlen = this->gyro_bias._encodeNoHash(buf, offset + pos, maxlen - pos);
    if(tlen < 0) return tlen; else pos += tlen;

    return pos;
}

__lcm_buffer_size imu_biases_t::_decodeNoHash(const void *buf, __lcm_buffer_size offset, __lcm_buffer_size maxlen)
{
    __lcm_buffer_size pos = 0, tlen;

    tlen = this->accel_bias._decodeNoHash(buf, offset + pos, maxlen - pos);
    if(tlen < 0) return tlen; else pos += tlen;

    tlen = this->gyro_bias._decodeNoHash(buf, offset + pos, maxlen - pos);
    if(tlen < 0) return tlen; else pos += tlen;

    return pos;
}

__lcm_buffer_size imu_biases_t::_getEncodedSizeNoHash() const
{
    __lcm_buffer_size enc_size = 0;
    enc_size += this->accel_bias._getEncodedSizeNoHash();
    enc_size += this->gyro_bias._getEncodedSizeNoHash();
    return enc_size;
}

bool imu_biases_t::operator==(const imu_biases_t& other) const {
  return (
          (accel_bias==other.accel_bias) && 
          (gyro_bias==other.gyro_bias));
}

bool imu_biases_t::operator!=(const imu_biases_t& other) const {
  return !(*this==other);
}

}

#endif