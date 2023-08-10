/*
    This file is part of Mitsuba, a physically based rendering system.

    Copyright (c) 2007-2014 by Wenzel Jakob and others.

    Mitsuba is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License Version 3
    as published by the Free Software Foundation.

    Mitsuba is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#include "pssmlt_sampler.h"
#include <iostream>
#include <fstream>

MTS_NAMESPACE_BEGIN

PSSMLTSampler::PSSMLTSampler(const PSSMLTConfiguration &config) : Sampler(Properties()) {
    m_random = new Random(config.seed);
    m_s1 = config.mutationSizeLow;
    m_s2 = config.mutationSizeHigh;
    m_width = config.width;
    m_height = config.height;
    m_sample_map_path = config.sampleMapPath;
    configure();
}

PSSMLTSampler::PSSMLTSampler(PSSMLTSampler *sampler) : Sampler(Properties()),
    m_random(sampler->m_random) {
    m_s1 = sampler->m_s1;
    m_s2 = sampler->m_s2;
    m_sample_map_path = sampler->m_sample_map_path;
    m_width = sampler->m_width;
    m_height = sampler->m_height;

    // Load sample map's CDF
    std::string sample_map_path_pdf = m_sample_map_path;
    if (sample_map_path_pdf.find("_pdf") != std::string::npos) {
        std::string sample_map_path_cdf = sample_map_path_pdf.replace(sample_map_path_pdf.find("_pdf"), 4, "_cdf");
        std::ifstream infile_cdf(sample_map_path_cdf.c_str(), std::ios::in | std::ios::binary);
        if (infile_cdf.is_open()) {
            adaptive = true;
            size_t size = m_width * m_height;
            sample_weight.resize(size);
            infile_cdf.read(reinterpret_cast<char*>(&sample_weight[0]), size * sizeof(double));
            infile_cdf.close();
        } 
    }
    
    configure();
}

PSSMLTSampler::PSSMLTSampler(Stream *stream, InstanceManager *manager)
    : Sampler(stream, manager) {
    m_random = static_cast<Random *>(manager->getInstance(stream));
    m_s1 = stream->readFloat();
    m_s2 = stream->readFloat();
    configure();
}

void PSSMLTSampler::serialize(Stream *stream, InstanceManager *manager) const {
    Sampler::serialize(stream, manager);
    manager->serialize(stream, m_random.get());
    stream->writeFloat(m_s1);
    stream->writeFloat(m_s2);
}

void PSSMLTSampler::configure() {
    m_logRatio = -math::fastlog(m_s2/m_s1);
    m_time = 0;
    m_largeStepTime = 0;
    m_largeStep = false;
    m_sampleIndex = 0;
    m_sampleCount = 0;
}

PSSMLTSampler::~PSSMLTSampler() { }

void PSSMLTSampler::accept() {
    if (m_largeStep)
        m_largeStepTime = m_time;
    m_time++;
    m_backup.clear();
    m_sampleIndex = 0;
}

void PSSMLTSampler::reset() {
    m_time = m_sampleIndex = m_largeStepTime = 0;
    m_u.clear();
}

void PSSMLTSampler::reject() {
    for (size_t i=0; i<m_backup.size(); ++i)
        m_u[m_backup[i].first] = m_backup[i].second;
    m_backup.clear();
    m_sampleIndex = 0;
}

Float PSSMLTSampler::primarySample(size_t i) {

    while (i >= m_u.size()) {
        m_u.push_back(SampleStruct(m_random->nextFloat()));
    }
    
    if (m_u[i].modify < m_time) {
        if (m_largeStep) {
            m_backup.push_back(std::pair<size_t, SampleStruct>(i, m_u[i]));
            m_u[i].modify = m_time;
            m_u[i].value = m_random->nextFloat();
            if (is_sensor) {
                current_adaptive = true;
            }
        } else {
            if (m_u[i].modify < m_largeStepTime) {
                m_u[i].modify = m_largeStepTime;
                m_u[i].value = m_random->nextFloat();
                if (is_sensor) {
                    current_adaptive = true;
                }
            }

            while (m_u[i].modify + 1 < m_time) {
                m_u[i].value = mutate(m_u[i].value);
                m_u[i].modify++;
            }

            m_backup.push_back(std::pair<size_t, SampleStruct>(i, m_u[i]));

            m_u[i].value = mutate(m_u[i].value);
            m_u[i].modify++;
        }
    }

    return m_u[i].value;
}

ref<Sampler> PSSMLTSampler::clone() {
    ref<PSSMLTSampler> sampler = new PSSMLTSampler(this);
    sampler->m_sampleCount = m_sampleCount;
    sampler->m_sampleIndex = m_sampleIndex;
    sampler->m_random = new Random(m_random);
    return sampler.get();
}

Float PSSMLTSampler::next1D() {
    return primarySample(m_sampleIndex++);
}

Point2 PSSMLTSampler::next2D() {

    /// Enforce a specific order of evaluation
    Float value1 = primarySample(m_sampleIndex++);
    Float value2 = primarySample(m_sampleIndex++);

    if (adaptive && current_adaptive) {
        current_adaptive = false;    
        Float target = m_random->nextFloat();
        int low = 0, high = m_width * m_height;
        while (low < high) {
            int mid = (low + high) >> 1;
            if (sample_weight[mid] < target)
                low = mid + 1;
            else
                high = mid;
        }
        Float row = low / m_width + m_random->nextFloat();
        Float col = low % m_width + m_random->nextFloat();
        Float row_norm = row / m_height;
        Float col_norm = col / m_width;
        m_u[m_sampleIndex - 2].value = col_norm;
        m_u[m_sampleIndex - 1].value = row_norm;
        return Point2(col_norm, row_norm);
    }

    return Point2(value1, value2);
}

std::string PSSMLTSampler::toString() const {
    std::ostringstream oss;
    oss << "PSSMLTSampler[" << endl
        << "  sampleCount = " << m_sampleCount << endl
        << "]";
    return oss.str();
}

MTS_IMPLEMENT_CLASS_S(PSSMLTSampler, false, Sampler)
MTS_NAMESPACE_END
