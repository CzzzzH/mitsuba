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

#include <mitsuba/render/bsdf.h>
#include <mitsuba/render/texture.h>
#include <mitsuba/hw/basicshader.h>
#include <mitsuba/core/warp.h>

MTS_NAMESPACE_BEGIN

class OpenRoomsMicrofacet : public BSDF {
public:
    OpenRoomsMicrofacet(const Properties &props)
        : BSDF(props) {
        /* For better compatibility with other models, support both
           'reflectance' and 'diffuseReflectance' as parameter names */
        m_albedo = new ConstantSpectrumTexture(props.getSpectrum("albedo", Spectrum(.5f)));
        m_roughness = new ConstantSpectrumTexture(props.getSpectrum("roughness", Spectrum(.0f)));
        m_normal = new ConstantSpectrumTexture(props.getSpectrum("normal", Spectrum(.0f)));
    }

    OpenRoomsMicrofacet(Stream *stream, InstanceManager *manager)
        : BSDF(stream, manager) {
        m_albedo = static_cast<Texture *>(manager->getInstance(stream));
        m_roughness = static_cast<Texture *>(manager->getInstance(stream));
        m_normal = static_cast<Texture *>(manager->getInstance(stream));
        configure();
    }

    void configure() {
        m_components.clear();
		m_components.push_back(EGlossyReflection | EFrontSide );
		m_components.push_back(EDiffuseReflection | EFrontSide );
		m_usesRayDifferentials  = false;
		BSDF::configure();
    }

    Spectrum eval(const BSDFSamplingRecord &bRec, EMeasure measure) const {
        
        Normal N;
        Frame frame = BSDF::getFrame(bRec.its);
        // m_normal->eval(bRec.its, false).toLinearRGB(N.x, N.y, N.z);
        // for (int i = 0; i < 3; ++i)
        //     N[i] = 2 * N[i] - 1;
        // N = normalize(frame.toWorld(N));
        N = frame.n;

        Vector V = bRec.wi;
        Vector L = bRec.wo;
        Vector H = normalize((L + V) / 2.0f);

        Float rough = m_roughness->eval(bRec.its)[0];
        Float alpha = rough * rough;
        Float k = (alpha + 2 * rough + 1) / 8.0;
        Float alpha2 = alpha * alpha;
        
        Float NoL = fmaxf(dot(N, L), 0);
        Float NoV = fmaxf(dot(N, V), 0);
        Float NoH = fmaxf(dot(N, H), 0);
        Float VoH = fmaxf(dot(V, H), 0);

        Float FMi = (-5.55473 * VoH - 6.98316) * VoH;
        Float F0 = 0.05;
        Spectrum frac0 = Spectrum(F0 + (1 - F0) * pow(2.0f, FMi));
        Spectrum frac = frac0 * alpha2;
        Float nom0 = NoH * NoH * (alpha2 - 1) + 1;
        Float nom1 = NoV * (1 - k) + k;
        Float nom2 = NoL * (1 - k) + k;
        Float nom = fmaxf(4 * M_PI * nom0 * nom0 * nom1 * nom2, 1e-14);
        Spectrum spec = frac / nom;
            
        Spectrum radiance = (m_albedo->eval(bRec.its) / M_PI + spec) * NoL;
        return radiance;
    }

    Float pdf(const BSDFSamplingRecord &bRec, EMeasure measure) const {
        
        Normal N;
        Frame frame = BSDF::getFrame(bRec.its);
        // m_normal->eval(bRec.its, false).toLinearRGB(N.x, N.y, N.z);
        // for (int i = 0; i < 3; ++i)
        //     N[i] = 2 * N[i] - 1;
        // N = normalize(frame.toWorld(N));
        N = frame.n;

        Vector V = bRec.wi;
        Vector L = bRec.wo;
        Float R = m_roughness->eval(bRec.its)[0];
        Float a2 = R * R * R * R;
        Vector H = normalize((L + V) / 2.0 );
        Float NoL = fmaxf(dot(N, L), 0);
        Float NoH = fmaxf(dot(N, H), 0);
        Float VoH = fmaxf(dot(V, H), 0);

        Float pdfLambertian = fmaxf(NoL / M_PI, 1e-14f);
        Float pdfSpecular = fmaxf((a2 * NoH) / fmaxf((4 * M_PI * (1 + (a2 - 1) * NoH) * (1 + (a2 - 1) * NoH) * VoH), 1e-14f), 1e-14f);
        return pdfLambertian * 0.5 + pdfSpecular * 0.5;
    }

    Spectrum sample(BSDFSamplingRecord &bRec, Float &pdf, const Point2 &sample) const {
        
        Point2 new_sample(sample);
        Normal N;
        Frame frame = BSDF::getFrame(bRec.its);
        // m_normal->eval(bRec.its, false).toLinearRGB(N.x, N.y, N.z);
        // for (int i = 0; i < 3; ++i)
        //     N[i] = 2 * N[i] - 1;
        // N = normalize(frame.toWorld(N));
        N = frame.n;

        Vector V = bRec.wi;
        Vector L;
        Float rough = m_roughness->eval(bRec.its)[0];
        
        Float alpha = rough * rough;
        Float k = (alpha + 2 * rough + 1) / 8.0;
        Float alpha2 = alpha * alpha;
    
        bool choseSpecular = true;
        if (new_sample.y < 0.5) {
            new_sample.y *= 2;
        } else {
            new_sample.y = (new_sample.y - 0.5) * 2;
            choseSpecular = false;
        }

        if (!choseSpecular) {
            bRec.sampledComponent = 1;
            bRec.sampledType = EDiffuseReflection;
            bRec.wo = warp::squareToCosineHemisphere(new_sample);
        } else {
            bRec.sampledComponent = 0;
            bRec.sampledType = EGlossyReflection;
            Float phi = 2 * M_PI * new_sample.x;
            Float cosTheta = sqrt((1 - new_sample.y) / (1 + (alpha2 - 1) * new_sample.y));
            Float sinTheta = sqrt(1 - cosTheta * cosTheta);
            Vector H = Vector(sinTheta * cos(phi), sinTheta * sin(phi), cosTheta);
            L = 2 * dot(V, H) * H - V;
            bRec.wo = L;
        }

        /* Guard against numerical imprecisions */
        pdf = OpenRoomsMicrofacet::pdf(bRec, ESolidAngle);

        if (pdf == 0)
            return Spectrum(0.0f);
        else
            return eval(bRec, ESolidAngle) / pdf;
    }

    Spectrum sample(BSDFSamplingRecord &bRec, const Point2 &sample) const {
        Float pdf;
        return OpenRoomsMicrofacet::sample(bRec, pdf, sample);
    }

    void addChild(const std::string &name, ConfigurableObject *child) {
        if (child->getClass()->derivesFrom(MTS_CLASS(Texture))) {
            if (name == "albedo")
                m_albedo = static_cast<Texture *>(child);
            else if (name == "roughness")
                m_roughness = static_cast<Texture *>(child);
            else if (name == "normal")
                m_normal = static_cast<Texture *>(child);
            else 
                BSDF::addChild(name, child);
        } else {
            BSDF::addChild(name, child);
        }
    }

    void serialize(Stream *stream, InstanceManager *manager) const {
        BSDF::serialize(stream, manager);
        manager->serialize(stream, m_albedo.get());
        manager->serialize(stream, m_roughness.get());
        manager->serialize(stream, m_normal.get());
    }

    Spectrum getDiffuseReflectance(const Intersection &its) const {
        return m_albedo->eval(its);
    }

    Float getRoughness(const Intersection &its, int component) const {
        return m_roughness->eval(its)[0];
    }

    std::string toString() const {
        std::ostringstream oss;
        oss << "OpenRoomsMicrofacet[" << endl
            << "  id = \"" << getID() << "\"," << endl
            << "  albedo = " << indent(m_albedo->toString()) << endl
            << "  roughness = " << indent(m_roughness->toString()) << endl
            << "  normal = " << indent(m_normal->toString()) << endl
            << "]";
        return oss.str();
    }

    Shader *createShader(Renderer *renderer) const;

    MTS_DECLARE_CLASS()
private:
    ref<Texture> m_albedo, m_normal, m_roughness;
};

// ================ Hardware shader implementation ================

class OpenRoomsMicrofacetShader : public Shader {
public:
	OpenRoomsMicrofacetShader(Renderer *renderer) :
		Shader(renderer, EBSDFShader) {
		m_flags = ETransparent;
	}

	void generateCode(std::ostringstream &oss,
			const std::string &evalName,
			const std::vector<std::string> &depNames) const {
		oss << "vec3 " << evalName << "(vec2 uv, vec3 wi, vec3 wo) {" << endl
			<< "    return vec3(0.0);" << endl
			<< "}" << endl;
		oss << "vec3 " << evalName << "_diffuse(vec2 uv, vec3 wi, vec3 wo) {" << endl
			<< "    return vec3(0.0);" << endl
			<< "}" << endl;
	}
	MTS_DECLARE_CLASS()
};

Shader *OpenRoomsMicrofacet::createShader(Renderer *renderer) const {
    return new OpenRoomsMicrofacetShader(renderer);
}

MTS_IMPLEMENT_CLASS(OpenRoomsMicrofacetShader, false, Shader)
MTS_IMPLEMENT_CLASS_S(OpenRoomsMicrofacet, false, BSDF)
MTS_EXPORT_PLUGIN(OpenRoomsMicrofacet, "OpenRooms microfacet BRDF")
MTS_NAMESPACE_END
