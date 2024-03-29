/** THIS IS AN AUTOMATICALLY GENERATED FILE.  DO NOT MODIFY
 * BY HAND!!
 *
 * Generated by lcm-gen
 *
 * From Source File: /Users/ryan/dev/symforce-zig-build/lcmtypes/symforce.lcm
 **/

#include <lcm/lcm_coretypes.h>

#ifndef __sym_optimization_stats_t_hpp__
#define __sym_optimization_stats_t_hpp__

#include <ostream>
#include <vector>
#include "lcmtypes/eigen_lcm/VectorXi.hpp"
#include "lcmtypes/sym/optimization_iteration_t.hpp"
#include "lcmtypes/sym/optimization_status_t.hpp"
#include "lcmtypes/sym/sparse_matrix_structure_t.hpp"

namespace sym
{

/// Debug stats for a full optimization run
class optimization_stats_t
{
    public:
        std::vector< ::sym::optimization_iteration_t > iterations;

        // Index into iterations of the best iteration (containing the optimal Values)
        int32_t best_index;

        // What was the result of the optimization? (did it converge, fail, etc.)
        ::sym::optimization_status_t status;

        // If status == FAILED, why?  This should be cast to the Optimizer::FailureReason enum
        // for the nonlinear solver you used.
        int32_t failure_reason;

        // The sparsity pattern of the jacobian, filled out if debug_stats=true
        ::sym::sparse_matrix_structure_t jacobian_sparsity;

        // The ordering used for the linear solver, filled out if debug_stats=true
        ::eigen_lcm::VectorXi linear_solver_ordering;

        // The sparsity pattern of the cholesky factor L, filled out if debug_stats=true
        ::sym::sparse_matrix_structure_t cholesky_factor_sparsity;

    public:
        optimization_stats_t() = default;

        /**
         * Member constructor
         */
        inline optimization_stats_t(
            const std::vector< ::sym::optimization_iteration_t >& iterations_arg,
            const int32_t& best_index_arg,
            const ::sym::optimization_status_t& status_arg,
            const int32_t& failure_reason_arg,
            const ::sym::sparse_matrix_structure_t& jacobian_sparsity_arg,
            const ::eigen_lcm::VectorXi& linear_solver_ordering_arg,
            const ::sym::sparse_matrix_structure_t& cholesky_factor_sparsity_arg
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

        using type_name_array_t = const char[21];

        inline static constexpr type_name_array_t* getTypeNameArrayPtr();

        /**
         * Returns "optimization_stats_t"
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
                if(fp->v == optimization_stats_t::getHash)
                    return 0;
            const __lcm_hash_ptr cp = { p, optimization_stats_t::getHash };

            uint64_t hash = 0xf934e9f17e183800LL +
                ::sym::optimization_iteration_t::_computeHash(&cp) +
         ::sym::optimization_status_t::_computeHash(&cp) +
         ::sym::sparse_matrix_structure_t::_computeHash(&cp) +
         ::eigen_lcm::VectorXi::_computeHash(&cp) +
         ::sym::sparse_matrix_structure_t::_computeHash(&cp);

            return (hash<<1) + ((hash>>63)&1);
        }

        // Comparison operators.
        inline bool operator==(const optimization_stats_t& other) const;
        inline bool operator!=(const optimization_stats_t& other) const;

        // Ability to print to standard streams as well as the fmt library.
        friend std::ostream& operator<<(std::ostream& stream, const optimization_stats_t& obj) {
#if defined(SKYMARSHAL_PRINTING_ENABLED)
            stream << "optimization_stats_t(";
            stream << "iterations=[";
            for (size_t i = 0; i < obj.iterations.size(); ++i) {
                stream << obj.iterations[i];
                if (i + 1 < obj.iterations.size()) {
                    stream << ", ";
                }
            }
            stream << "]" << ", ";
            stream << "best_index=" << obj.best_index << ", ";
            stream << "status=" << obj.status << ", ";
            stream << "failure_reason=" << obj.failure_reason << ", ";
            stream << "jacobian_sparsity=" << obj.jacobian_sparsity << ", ";
            stream << "linear_solver_ordering=" << obj.linear_solver_ordering << ", ";
            stream << "cholesky_factor_sparsity=" << obj.cholesky_factor_sparsity;
            stream << ")";
#else
            stream << "<FORMATTING DISABLED>";
#endif
            return stream;
        }
};

optimization_stats_t::optimization_stats_t(
    const std::vector< ::sym::optimization_iteration_t >& iterations_arg,
    const int32_t& best_index_arg,
    const ::sym::optimization_status_t& status_arg,
    const int32_t& failure_reason_arg,
    const ::sym::sparse_matrix_structure_t& jacobian_sparsity_arg,
    const ::eigen_lcm::VectorXi& linear_solver_ordering_arg,
    const ::sym::sparse_matrix_structure_t& cholesky_factor_sparsity_arg
) : iterations(iterations_arg),
    best_index(best_index_arg),
    status(status_arg),
    failure_reason(failure_reason_arg),
    jacobian_sparsity(jacobian_sparsity_arg),
    linear_solver_ordering(linear_solver_ordering_arg),
    cholesky_factor_sparsity(cholesky_factor_sparsity_arg) {}

__lcm_buffer_size optimization_stats_t::encode(void *buf, __lcm_buffer_size offset, __lcm_buffer_size maxlen) const
{
    __lcm_buffer_size pos = 0, tlen;
    int64_t hash = (int64_t)getHash();

    tlen = __int64_t_encode_array(buf, offset + pos, maxlen - pos, &hash, 1);
    if(tlen < 0) return tlen; else pos += tlen;

    tlen = this->_encodeNoHash(buf, offset + pos, maxlen - pos);
    if (tlen < 0) return tlen; else pos += tlen;

    return pos;
}

__lcm_buffer_size optimization_stats_t::decode(const void *buf, __lcm_buffer_size offset, __lcm_buffer_size maxlen)
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

__lcm_buffer_size optimization_stats_t::getEncodedSize() const
{
    return 8 + _getEncodedSizeNoHash();
}

int64_t optimization_stats_t::getHash()
{
    static int64_t hash = _computeHash(NULL);
    return hash;
}

constexpr optimization_stats_t::type_name_array_t* optimization_stats_t::getTypeNameArrayPtr() {
    return &"optimization_stats_t";
}

constexpr const char* optimization_stats_t::getTypeName()
{
    return *optimization_stats_t::getTypeNameArrayPtr();
}

constexpr const char * optimization_stats_t::getPackageName()
{
    return "sym";
}

__lcm_buffer_size optimization_stats_t::_encodeNoHash(void *buf, __lcm_buffer_size offset, __lcm_buffer_size maxlen) const
{
    __lcm_buffer_size pos = 0, tlen;

    int32_t v_num_iterations = this->iterations.size();
    tlen = __int32_t_encode_array(buf, offset + pos, maxlen - pos, &v_num_iterations, 1);
    if(tlen < 0) return tlen; else pos += tlen;

    for (__lcm_buffer_size a0 = 0; a0 < v_num_iterations; a0++) {
        tlen = this->iterations[a0]._encodeNoHash(buf, offset + pos, maxlen - pos);
        if(tlen < 0) return tlen; else pos += tlen;
    }

    tlen = __int32_t_encode_array(buf, offset + pos, maxlen - pos, &this->best_index, 1);
    if(tlen < 0) return tlen; else pos += tlen;

    tlen = this->status._encodeNoHash(buf, offset + pos, maxlen - pos);
    if(tlen < 0) return tlen; else pos += tlen;

    tlen = __int32_t_encode_array(buf, offset + pos, maxlen - pos, &this->failure_reason, 1);
    if(tlen < 0) return tlen; else pos += tlen;

    tlen = this->jacobian_sparsity._encodeNoHash(buf, offset + pos, maxlen - pos);
    if(tlen < 0) return tlen; else pos += tlen;

    tlen = this->linear_solver_ordering._encodeNoHash(buf, offset + pos, maxlen - pos);
    if(tlen < 0) return tlen; else pos += tlen;

    tlen = this->cholesky_factor_sparsity._encodeNoHash(buf, offset + pos, maxlen - pos);
    if(tlen < 0) return tlen; else pos += tlen;

    return pos;
}

__lcm_buffer_size optimization_stats_t::_decodeNoHash(const void *buf, __lcm_buffer_size offset, __lcm_buffer_size maxlen)
{
    __lcm_buffer_size pos = 0, tlen;

    int32_t v_num_iterations;
    tlen = __int32_t_decode_array(buf, offset + pos, maxlen - pos, &v_num_iterations, 1);
    if(tlen < 0) return tlen; else pos += tlen;

    this->iterations.resize(v_num_iterations);
    for (__lcm_buffer_size a0 = 0; a0 < v_num_iterations; a0++) {
        tlen = this->iterations[a0]._decodeNoHash(buf, offset + pos, maxlen - pos);
        if(tlen < 0) return tlen; else pos += tlen;
    }

    tlen = __int32_t_decode_array(buf, offset + pos, maxlen - pos, &this->best_index, 1);
    if(tlen < 0) return tlen; else pos += tlen;

    tlen = this->status._decodeNoHash(buf, offset + pos, maxlen - pos);
    if(tlen < 0) return tlen; else pos += tlen;

    tlen = __int32_t_decode_array(buf, offset + pos, maxlen - pos, &this->failure_reason, 1);
    if(tlen < 0) return tlen; else pos += tlen;

    tlen = this->jacobian_sparsity._decodeNoHash(buf, offset + pos, maxlen - pos);
    if(tlen < 0) return tlen; else pos += tlen;

    tlen = this->linear_solver_ordering._decodeNoHash(buf, offset + pos, maxlen - pos);
    if(tlen < 0) return tlen; else pos += tlen;

    tlen = this->cholesky_factor_sparsity._decodeNoHash(buf, offset + pos, maxlen - pos);
    if(tlen < 0) return tlen; else pos += tlen;

    return pos;
}

__lcm_buffer_size optimization_stats_t::_getEncodedSizeNoHash() const
{
    __lcm_buffer_size enc_size = 0;
    enc_size += __int32_t_encoded_array_size(NULL, 1);
    for (__lcm_buffer_size a0 = 0; a0 < this->iterations.size(); a0++) {
        enc_size += this->iterations[a0]._getEncodedSizeNoHash();
    }
    enc_size += __int32_t_encoded_array_size(NULL, 1);
    enc_size += this->status._getEncodedSizeNoHash();
    enc_size += __int32_t_encoded_array_size(NULL, 1);
    enc_size += this->jacobian_sparsity._getEncodedSizeNoHash();
    enc_size += this->linear_solver_ordering._getEncodedSizeNoHash();
    enc_size += this->cholesky_factor_sparsity._getEncodedSizeNoHash();
    return enc_size;
}

bool optimization_stats_t::operator==(const optimization_stats_t& other) const {
  return (
          (iterations==other.iterations) && 
          (best_index==other.best_index) && 
          (status==other.status) && 
          (failure_reason==other.failure_reason) && 
          (jacobian_sparsity==other.jacobian_sparsity) && 
          (linear_solver_ordering==other.linear_solver_ordering) && 
          (cholesky_factor_sparsity==other.cholesky_factor_sparsity));
}

bool optimization_stats_t::operator!=(const optimization_stats_t& other) const {
  return !(*this==other);
}

}

#endif