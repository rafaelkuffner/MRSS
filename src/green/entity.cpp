
#include "entity.hpp"

#include <atomic>

green::Entity::Entity() {
	static std::atomic<int> nextid{1};
	m_id = nextid++;
}
