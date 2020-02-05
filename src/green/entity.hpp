
#pragma once

#ifndef GREEN_ENTITY_HPP
#define GREEN_ENTITY_HPP

#include <cgu/util.hpp>

namespace green {

	class Entity {
	private:
		Entity(const Entity &) = delete;
		Entity & operator=(const Entity &) = delete;

	public:
		Entity() {}

		virtual bool dead() const = 0;

		// model to world
		virtual glm::mat4 transform() const = 0;

		virtual void draw(const glm::mat4 &view, const glm::mat4 &proj, float zfar) = 0;

		virtual ~Entity() {}
	};

}

#endif
