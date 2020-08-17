
#pragma once

#ifndef GREEN_ENTITY_HPP
#define GREEN_ENTITY_HPP

#include <cgu/util.hpp>

#include "saliency.hpp"

namespace green {

	struct entity_selection {
		int hover_entity = -1;
		int hover_vertex = -1;
		int select_entity = -1;
		int select_vertex = -1;
		float hover_entity_dist = 0;
		float hover_vertex_dist = 0;
	};

	class Entity {
	private:
		Entity(const Entity &) = delete;
		Entity & operator=(const Entity &) = delete;

		int m_id;

	public:
		Entity();

		int id() const {
			return m_id;
		}

		virtual std::string name() const = 0;

		virtual void try_kill() = 0;

		virtual bool dead() const = 0;

		virtual void move_by(const glm::vec3 &d) {}

		// model to world
		virtual glm::mat4 transform() const = 0;

		virtual void draw(const glm::mat4 &view, const glm::mat4 &proj, float zfar, bool draw_scene) = 0;

		virtual std::future<saliency_result> compute_saliency_async(const saliency_user_params &uparams, saliency_progress &progress) {
			return {};
		}

		virtual ~Entity() {}
	};

}

#endif
